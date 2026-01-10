<#
.SYNOPSIS
    Parallel Simulation Runner for Chencode Simulator
    
.DESCRIPTION
    Automatically detects CPU cores and runs multiple simulation instances
    in parallel, then merges the results into a single CSV file.
    
.PARAMETER Decoder
    Decoder type: 0=Uncoded, 1=HardViterbi, 2=SoftViterbi, 3=BCJR, 4=Turbo
    
.PARAMETER StartSNR
    Starting Eb/N0 in dB

.PARAMETER EndSNR
    Ending Eb/N0 in dB

.PARAMETER Step
    SNR step size in dB

.PARAMETER Frames
    Number of frames per SNR point

.PARAMETER Workers
    Number of parallel workers (default: auto-detect based on CPU cores)
    
.EXAMPLE
    .\parallel_sim.ps1 -Decoder 4 -StartSNR -1.0 -EndSNR 4.0 -Step 0.5 -Frames 50000
#>

param(
    [Parameter(Mandatory=$true)]
    [ValidateRange(0, 4)]
    [int]$Decoder,
    
    [Parameter(Mandatory=$false)]
    [ValidateRange(0, 1)]
    [int]$TurboType = 0,  # 0=(7,5)_8, 1=CCSDS
    
    [Parameter(Mandatory=$true)]
    [float]$StartSNR,
    
    [Parameter(Mandatory=$true)]
    [float]$EndSNR,
    
    [Parameter(Mandatory=$true)]
    [float]$Step,
    
    [Parameter(Mandatory=$true)]
    [long]$Frames,
    
    [int]$Workers = 0
)

# Configuration
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$ProjectRoot = Split-Path -Parent $ScriptDir
$ExePath = Join-Path $ProjectRoot "chencode_sim.exe"
$OutputDir = Join-Path $ProjectRoot "output"
$Timestamp = Get-Date -Format "yyyyMMdd_HHmmss"

$DecoderNames = @(
    "Uncoded_BPSK",
    "Hard_Viterbi",
    "Soft_Viterbi",
    "BCJR_MAP",
    "Turbo_LogMAP"
)

# System Detection
function Get-SystemInfo {
    $cpu = Get-CimInstance -ClassName Win32_Processor
    $cpuName = $cpu.Name.Trim()
    $physicalCores = $cpu.NumberOfCores
    $logicalCores = $cpu.NumberOfLogicalProcessors
    $memory = [math]::Round((Get-CimInstance -ClassName Win32_ComputerSystem).TotalPhysicalMemory / 1GB, 1)
    
    return @{
        CpuName = $cpuName
        PhysicalCores = $physicalCores
        LogicalCores = $logicalCores
        MemoryGB = $memory
    }
}

# SNR Point Distribution
function Get-SNRPoints {
    param([float]$Start, [float]$End, [float]$StepSize)
    
    $points = @()
    $snr = $Start
    while ($snr -le $End + 0.001) {
        $points += [math]::Round($snr, 2)
        $snr += $StepSize
    }
    return $points
}

function Split-SNRPoints {
    param(
        [float[]]$Points,
        [int]$NumWorkers
    )
    
    $chunks = @()
    $chunkSize = [math]::Ceiling($Points.Count / $NumWorkers)
    
    for ($i = 0; $i -lt $NumWorkers; $i++) {
        $startIdx = $i * $chunkSize
        $endIdx = [math]::Min(($i + 1) * $chunkSize - 1, $Points.Count - 1)
        
        if ($startIdx -le $endIdx) {
            $chunks += ,@{
                WorkerId = $i
                Points = $Points[$startIdx..$endIdx]
                StartSNR = $Points[$startIdx]
                EndSNR = $Points[$endIdx]
            }
        }
    }
    return $chunks
}

# Progress Monitoring
function Get-PartialProgress {
    param([string]$CsvFile)
    
    if (-not (Test-Path $CsvFile)) {
        return 0
    }
    
    try {
        $content = Get-Content $CsvFile -ErrorAction SilentlyContinue
        $dataLines = ($content | Where-Object { $_ -notmatch "^#" -and $_.Trim() -ne "" -and $_ -notmatch "Eb_N0" }).Count
        return $dataLines
    } catch {
        return 0
    }
}

function Show-Progress {
    param(
        [System.Diagnostics.Process[]]$Processes,
        [hashtable[]]$Chunks,
        [string[]]$OutputFiles
    )
    
    $allDone = $false
    $startTime = Get-Date
    
    while (-not $allDone) {
        Start-Sleep -Seconds 3
        
        $totalExpected = 0
        $totalCompleted = 0
        
        Write-Host ""
        $elapsed = (Get-Date) - $startTime
        $elapsedStr = "{0:hh\:mm\:ss}" -f $elapsed
        Write-Host "  Progress (Elapsed: $elapsedStr):" -ForegroundColor Cyan
        
        for ($i = 0; $i -lt $Chunks.Count; $i++) {
            $expected = $Chunks[$i].Points.Count
            $completed = Get-PartialProgress -CsvFile $OutputFiles[$i]
            $totalExpected += $expected
            $totalCompleted += $completed
            
            $proc = $Processes[$i]
            if ($proc.HasExited) {
                $status = "Done"
            } else {
                $status = "Running"
            }
            
            $pct = 0
            if ($expected -gt 0) {
                $pct = [math]::Round(100 * $completed / $expected)
            }
            
            $snrS = $Chunks[$i].StartSNR
            $snrE = $Chunks[$i].EndSNR
            Write-Host "    Worker $i : SNR $snrS to $snrE : $completed / $expected pts ($pct pct) - $status"
        }
        
        $overallPct = 0
        if ($totalExpected -gt 0) {
            $overallPct = [math]::Round(100 * $totalCompleted / $totalExpected)
        }
        Write-Host "    TOTAL: $totalCompleted / $totalExpected pts ($overallPct pct)" -ForegroundColor Yellow
        
        $running = ($Processes | Where-Object { -not $_.HasExited }).Count
        $allDone = ($running -eq 0)
    }
    
    return $elapsed
}

# CSV Merging
function Merge-CSVFiles {
    param(
        [string[]]$InputFiles,
        [string]$OutputFile,
        [string]$DecoderName,
        [long]$TotalFrames
    )
    
    $allData = @()
    
    foreach ($file in $InputFiles) {
        if (Test-Path $file) {
            $content = Get-Content $file
            foreach ($line in $content) {
                # Skip comments, empty lines, and header rows
                if ($line -notmatch "^#" -and $line.Trim() -ne "" -and $line -notmatch "Eb_N0" -and $line -notmatch "SNR_dB") {
                    $allData += $line
                }
            }
        }
    }
    
    $sortedData = $allData | Sort-Object { [float]($_ -split ",")[0] }
    
    $header = @(
        "# Channel Coding Simulation Results (Merged)",
        "# Decoder: $DecoderName",
        "# Total Frames per SNR: $TotalFrames",
        "# Generated: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')",
        "#",
        "Eb_N0,Bit_Errors,Total_Bits,BER,Frame_Errors,Total_Frames,FER"
    )
    
    $header | Out-File -FilePath $OutputFile -Encoding UTF8
    $sortedData | Out-File -FilePath $OutputFile -Encoding UTF8 -Append
    
    return $sortedData.Count
}

# === Main Execution ===

Write-Host ""
Write-Host "  ================================================================" -ForegroundColor Cyan
Write-Host "           CHENCODE Parallel Simulation Runner" -ForegroundColor Cyan
Write-Host "  ================================================================" -ForegroundColor Cyan

# Check executable exists
if (-not (Test-Path $ExePath)) {
    Write-Host "  [Error] Executable not found: $ExePath" -ForegroundColor Red
    Write-Host "  Please run build_msvc.bat first." -ForegroundColor Red
    exit 1
}

# Get system info
$sysInfo = Get-SystemInfo
Write-Host ""
Write-Host "  [System] $($sysInfo.CpuName)" -ForegroundColor Green
Write-Host "           $($sysInfo.PhysicalCores) cores, $($sysInfo.LogicalCores) threads, $($sysInfo.MemoryGB) GB RAM"

# Determine worker count
if ($Workers -eq 0) {
    $Workers = $sysInfo.PhysicalCores
}
Write-Host "  [Config] Using $Workers parallel workers" -ForegroundColor Green

# Calculate SNR points
$snrPoints = Get-SNRPoints -Start $StartSNR -End $EndSNR -StepSize $Step
Write-Host "  [Config] $($snrPoints.Count) SNR points: $StartSNR to $EndSNR dB (step $Step)"
Write-Host "  [Config] $Frames frames per SNR point"
Write-Host "  [Config] Decoder: $($DecoderNames[$Decoder])"

# Ensure output directory exists
if (-not (Test-Path $OutputDir)) {
    New-Item -ItemType Directory -Path $OutputDir -Force | Out-Null
}

# Split work
$chunks = Split-SNRPoints -Points $snrPoints -NumWorkers $Workers
$actualWorkers = $chunks.Count

Write-Host ""
Write-Host "  [Info] Work distribution:" -ForegroundColor Yellow
for ($w = 0; $w -lt $chunks.Count; $w++) {
    $c = $chunks[$w]
    $ss = $c.StartSNR
    $se = $c.EndSNR
    $np = $c.Points.Count
    Write-Host "    Worker $w : SNR $ss to $se ($np points)"
}

# Prepare output files and processes
$outputFiles = @()
$processes = @()

Write-Host ""
Write-Host "  [Start] Launching $actualWorkers workers..." -ForegroundColor Cyan

for ($i = 0; $i -lt $actualWorkers; $i++) {
    $chunk = $chunks[$i]
    $outFile = Join-Path $OutputDir "part_${i}_${Timestamp}.csv"
    $outputFiles += $outFile
    
    $seed = 10000 + $i * 1000 + (Get-Random -Maximum 999)
    
    # Build command arguments (include --turbo-type for Turbo decoder)
    $turboArg = ""
    if ($Decoder -eq 4) {
        $turboArg = " --turbo-type $TurboType"
    }
    $cmdArgs = "--batch --decoder $Decoder$turboArg --snr $($chunk.StartSNR) $($chunk.EndSNR) $Step --frames $Frames --output `"$outFile`" --seed $seed --quiet"
    
    $proc = Start-Process -FilePath $ExePath -ArgumentList $cmdArgs -PassThru -WindowStyle Hidden
    $processes += $proc
    
    Write-Host "    Worker $i : PID $($proc.Id) started"
}

# Monitor progress
$elapsed = Show-Progress -Processes $processes -Chunks $chunks -OutputFiles $outputFiles

# Wait for all processes to complete
$processes | ForEach-Object { $_.WaitForExit() }

Write-Host ""
Write-Host "  [Done] All workers completed!" -ForegroundColor Green

# Merge results
$mergedFile = Join-Path $OutputDir ("ber_" + $DecoderNames[$Decoder] + "_" + $Timestamp + "_merged.csv")
Write-Host "  [Merge] Combining results..." -ForegroundColor Cyan

$dataPoints = Merge-CSVFiles -InputFiles $outputFiles -OutputFile $mergedFile -DecoderName $DecoderNames[$Decoder] -TotalFrames $Frames

# Clean up partial files
foreach ($file in $outputFiles) {
    if (Test-Path $file) {
        Remove-Item $file -Force
    }
}

# Summary
$elapsedStr = "{0:hh\:mm\:ss}" -f $elapsed
Write-Host ""
Write-Host "  ================================================================" -ForegroundColor Green
Write-Host "    Simulation Complete!" -ForegroundColor Green
Write-Host "  ================================================================" -ForegroundColor Green
Write-Host "    Total time:    $elapsedStr"
Write-Host "    Data points:   $dataPoints"
Write-Host "    Output file:   $mergedFile"
Write-Host ""

<#
.SYNOPSIS
    Parallel Simulation Runner for Chencode Simulator
    
.DESCRIPTION
    Automatically detects CPU cores and runs multiple simulation instances
    in parallel, then merges the results into a single CSV file.
    
.PARAMETER Decoder
    Decoder type: 0=Uncoded, 1=HardViterbi, 2=SoftViterbi, 3=BCJR, 4=Turbo
    
.PARAMETER StartSNR
    Starting Eb/N0 in dB (e.g., -1.0)
    
.PARAMETER EndSNR
    Ending Eb/N0 in dB (e.g., 4.0)
    
.PARAMETER Step
    SNR step size in dB (e.g., 0.5)
    
.PARAMETER Frames
    Number of frames per SNR point (e.g., 50000)
    
.PARAMETER Workers
    Number of parallel workers (default: auto-detect based on CPU cores)
    
.EXAMPLE
    .\parallel_sim.ps1 -Decoder 4 -StartSNR -1.0 -EndSNR 4.0 -Step 0.5 -Frames 50000
#>

param(
    [Parameter(Mandatory=$true)]
    [ValidateRange(0, 4)]
    [int]$Decoder,
    
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

# =================================================================
# Configuration
# =================================================================

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$ProjectRoot = Split-Path -Parent $ScriptDir
$ExePath = Join-Path $ProjectRoot "chencode_sim.exe"
$OutputDir = Join-Path $ProjectRoot "output"
$Timestamp = Get-Date -Format "yyyyMMdd_HHmmss"

# Decoder names for display
$DecoderNames = @(
    "Uncoded BPSK",
    "Hard Viterbi (CC R=1/2)",
    "Soft Viterbi (CC R=1/2)",
    "BCJR / MAP (CC R=1/2)",
    "Turbo (Log-MAP, PCCC R=1/3)"
)

# =================================================================
# System Detection
# =================================================================

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

# =================================================================
# SNR Point Distribution
# =================================================================

function Get-SNRPoints {
    param([float]$Start, [float]$End, [float]$Step)
    
    $points = @()
    $snr = $Start
    while ($snr -le $End + 0.001) {
        $points += [math]::Round($snr, 2)
        $snr += $Step
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

# =================================================================
# Progress Monitoring
# =================================================================

function Get-PartialProgress {
    param([string]$CsvFile, [int]$ExpectedPoints)
    
    if (-not (Test-Path $CsvFile)) {
        return 0
    }
    
    try {
        $lines = Get-Content $CsvFile -ErrorAction SilentlyContinue | Where-Object { $_ -notmatch "^#" -and $_.Trim() -ne "" }
        $dataLines = ($lines | Measure-Object).Count - 1  # Subtract header
        if ($dataLines -lt 0) { $dataLines = 0 }
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
        Start-Sleep -Seconds 2
        
        $totalExpected = 0
        $totalCompleted = 0
        $statusLines = @()
        
        for ($i = 0; $i -lt $Chunks.Count; $i++) {
            $expected = $Chunks[$i].Points.Count
            $completed = Get-PartialProgress -CsvFile $OutputFiles[$i] -ExpectedPoints $expected
            $totalExpected += $expected
            $totalCompleted += $completed
            
            $proc = $Processes[$i]
            $status = if ($proc.HasExited) { "Done" } else { "Running" }
            $pct = if ($expected -gt 0) { [math]::Round(100 * $completed / $expected) } else { 0 }
            
            $statusLines += "    [Worker $i] SNR $($Chunks[$i].StartSNR) to $($Chunks[$i].EndSNR): $completed/$expected points ($pct%) - $status"
        }
        
        # Clear previous lines and print status
        $elapsed = (Get-Date) - $startTime
        $elapsedStr = "{0:hh\:mm\:ss}" -f $elapsed
        
        Write-Host "`r`n  Progress (Elapsed: $elapsedStr):" -ForegroundColor Cyan
        foreach ($line in $statusLines) {
            Write-Host $line
        }
        
        $overallPct = if ($totalExpected -gt 0) { [math]::Round(100 * $totalCompleted / $totalExpected) } else { 0 }
        Write-Host "    [Total] $totalCompleted/$totalExpected points ($overallPct%)" -ForegroundColor Yellow
        
        # Check if all processes are done
        $allDone = ($Processes | Where-Object { -not $_.HasExited } | Measure-Object).Count -eq 0
        
        if (-not $allDone) {
            # Move cursor up to overwrite progress
            $linesToMove = $Chunks.Count + 3
            Write-Host "`e[$($linesToMove)A" -NoNewline
        }
    }
    
    return $elapsed
}

# =================================================================
# CSV Merging
# =================================================================

function Merge-CSVFiles {
    param(
        [string[]]$InputFiles,
        [string]$OutputFile,
        [string]$DecoderName,
        [int]$TotalFrames
    )
    
    $allData = @()
    
    foreach ($file in $InputFiles) {
        if (Test-Path $file) {
            $content = Get-Content $file
            foreach ($line in $content) {
                if ($line -notmatch "^#" -and $line.Trim() -ne "" -and $line -notmatch "Eb_N0") {
                    $allData += $line
                }
            }
        }
    }
    
    # Sort by SNR (first column)
    $sortedData = $allData | Sort-Object { [float]($_ -split ",")[0] }
    
    # Write merged file
    $header = @(
        "# Channel Coding Simulation Results (Merged)",
        "# Decoder: $DecoderName",
        "# Total Frames per SNR: $TotalFrames",
        "# Generated: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')",
        "#",
        "Eb_N0,Bit_Errors,Total_Bits,BER"
    )
    
    $header | Out-File -FilePath $OutputFile -Encoding UTF8
    $sortedData | Out-File -FilePath $OutputFile -Encoding UTF8 -Append
    
    return $sortedData.Count
}

# =================================================================
# Main Execution
# =================================================================

# Print CHENCODE ASCII art banner with gradient colors
Write-Host ""
Write-Host ""
Write-Host "   " -NoNewline; Write-Host "*" -ForegroundColor DarkYellow -NoNewline; Write-Host " Welcome to " -NoNewline; Write-Host "Chencode Parallel Runner" -ForegroundColor White
Write-Host ""

# CHENCODE ASCII Art (matching main.c style)
Write-Host "    ██████╗██╗  ██╗███████╗███╗   ██╗ ██████╗ ██████╗ ██████╗ ███████╗" -ForegroundColor DarkYellow
Write-Host "   ██╔════╝██║  ██║██╔════╝████╗  ██║██╔════╝██╔═══██╗██╔══██╗██╔════╝" -ForegroundColor Yellow
Write-Host "   ██║     ███████║█████╗  ██╔██╗ ██║██║     ██║   ██║██║  ██║█████╗  " -ForegroundColor Green
Write-Host "   ██║     ██╔══██║██╔══╝  ██║╚██╗██║██║     ██║   ██║██║  ██║██╔══╝  " -ForegroundColor Cyan
Write-Host "   ╚██████╗██║  ██║███████╗██║ ╚████║╚██████╗╚██████╔╝██████╔╝███████╗" -ForegroundColor Blue
Write-Host "    ╚═════╝╚═╝  ╚═╝╚══════╝╚═╝  ╚═══╝ ╚═════╝ ╚═════╝ ╚═════╝ ╚══════╝" -ForegroundColor DarkCyan

Write-Host ""
Write-Host "                    Channel Coding Parallel Simulator"
Write-Host "                 (7,5)_8 Convolutional & PCCC Turbo Codes"
Write-Host ""
Write-Host "  ================================================================" -ForegroundColor Cyan
Write-Host "              Multi-Process Parallel Execution Mode" -ForegroundColor Cyan
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
$snrPoints = Get-SNRPoints -Start $StartSNR -End $EndSNR -Step $Step
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
foreach ($chunk in $chunks) {
    Write-Host "    Worker $($chunk.WorkerId): SNR $($chunk.StartSNR) to $($chunk.EndSNR) ($($chunk.Points.Count) points)"
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
    
    $args = "--batch --decoder $Decoder --snr $($chunk.StartSNR) $($chunk.EndSNR) $Step --frames $Frames --output `"$outFile`" --seed $seed --quiet"
    
    $proc = Start-Process -FilePath $ExePath -ArgumentList $args -PassThru -WindowStyle Hidden
    $processes += $proc
    
    Write-Host "    [Worker $i] PID $($proc.Id) started"
}

# Monitor progress
Write-Host ""
$elapsed = Show-Progress -Processes $processes -Chunks $chunks -OutputFiles $outputFiles

# Wait for all processes to complete
$processes | ForEach-Object { $_.WaitForExit() }

Write-Host ""
Write-Host "  [Done] All workers completed!" -ForegroundColor Green

# Merge results
$mergedFile = Join-Path $OutputDir "ber_$($DecoderNames[$Decoder] -replace '[^a-zA-Z0-9]', '_')_${Timestamp}_merged.csv"
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
Write-Host "    MATLAB import:" -ForegroundColor Yellow
Write-Host "    data = readtable('$mergedFile', 'CommentStyle', '#');"
Write-Host "    semilogy(data.Eb_N0, data.BER, '-o');"
Write-Host ""

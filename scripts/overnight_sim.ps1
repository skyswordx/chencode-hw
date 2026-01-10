<#
.SYNOPSIS
    Overnight Batch Simulation Scheduler
    
.DESCRIPTION
    Runs multiple simulation configurations in a round-robin fashion to ensure even coverage.
    Supports time limits, automatic merging, and progress reporting.
    
.PARAMETER ConfigFile
    Path to the JSON configuration file defining tasks.
    
.PARAMETER MaxHours
    Maximum execution time in hours.
    
.PARAMETER BatchSize
    Number of frames to run per task per round.
    
.PARAMETER DryRun
    If set, only lists tasks without executing them.
#>

param(
    [string]$ConfigFile = "sim_tasks.json",
    [double]$MaxHours = 0, # Set to 0 for infinite loop
    [int]$BatchSize = 10000,
    [string]$OutputDirName = "output",
    [switch]$DryRun
)

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$ProjectRoot = Split-Path -Parent $ScriptDir
$OutputDir = Join-Path $ProjectRoot $OutputDirName
$ConfigPath = Join-Path $ScriptDir $ConfigFile
$ParallelSimPath = Join-Path $ScriptDir "parallel_sim.ps1"

# === 1. Setup & Validation ===

if (-not (Test-Path $ConfigPath)) {
    Write-Error "Configuration file not found: $ConfigPath"
    exit 1
}

# Ensure output directory exists
if (-not (Test-Path $OutputDir)) {
    New-Item -ItemType Directory -Path $OutputDir -Force | Out-Null
}

$config = Get-Content $ConfigPath | Out-String | ConvertFrom-Json
$tasks = $config.tasks
$startTime = Get-Date

Write-Host "==========================================================" -ForegroundColor Cyan
Write-Host "   Overnight Batch Simulation Scheduler" -ForegroundColor Cyan
Write-Host "==========================================================" -ForegroundColor Cyan
Write-Host "Start Time: $($startTime)"
if ($MaxHours -gt 0) {
    Write-Host "Max Duration: $MaxHours hours"
} else {
    Write-Host "Max Duration: Infinite (Stop with Ctrl+C)"
}
Write-Host "Tasks: $($tasks.Count)"
Write-Host "Batch Size: $BatchSize frames"
Write-Host "Config: $ConfigFile"
Write-Host ""

# State tracking
$taskState = @{}
foreach ($t in $tasks) {
    $taskState[$t.id] = @{
        CompletedCycles = 0
        CompletedFrames = 0 # Track frames too
        Status = "Pending"
        LastOutputFile = ""
        TotalDuration = [timespan]::Zero
    }
}

if ($DryRun) {
    Write-Host "[Dry Run] Task List:" -ForegroundColor Yellow
    $tasks | Format-Table id, name, ccsdsK, ccsdsRate, targetFrames
    exit 0
}

# === 2. Main Execution Loop ===

$timeLimit = [timespan]::FromHours($MaxHours)
$globalIndex = 0
$pendingTasks = $tasks # Initialize pendingTasks for round-robin

while ($true) {
    # Check global time limit
    if ($MaxHours -gt 0) {
        $elapsed = (Get-Date) - $startTime
        if ($elapsed -ge $timeLimit) {
            Write-Host "`n[Time Limit Reached] Stopping after $elapsed..." -ForegroundColor Yellow
            break
        }
    }
    
    # Check if all tasks are done (ONLY if MaxHours > 0)
    # If MaxHours is 0 (infinite), we ignore target frames and just keep going round-robin
    # Check if all tasks are done
    # For now, we assume infinite loop or MaxHours governs everything
    # If using strictly MaxHours=0, we never stop until manual break.
    
    # Select next task (Round Robin)
    $currentTask = $pendingTasks[$globalIndex % $pendingTasks.Count]
    $globalIndex++
    
    $state = $taskState[$currentTask.id]
    
    $state['Status'] = "Running"
    
    if ($MaxHours -gt 0 -and $currentTask.targetFrames -gt 0) {
        Write-Host "`n[Running] Task: $($currentTask.name) (Completed: $($state['CompletedFrames'])/$($currentTask.targetFrames) frames)" -ForegroundColor Cyan
        $progress = [math]::Round(($state['CompletedFrames'] / $currentTask.targetFrames) * 100, 1)
        Write-Host "  > Task Progress: $progress% complete." -ForegroundColor Green
    } else {
        Write-Host "`n[Running] Task: $($currentTask.name) (Cycle: $($state['CompletedCycles'] + 1))" -ForegroundColor Cyan
        Write-Host "  > Task Progress: Continuous mode (Cycle: $($state['CompletedCycles']))" -ForegroundColor Green
    }
    Write-Host "  > Running SNR sweep with $BatchSize frames per point..."
    
    $batchStart = Get-Date
    
    $simArgs = @{
        Decoder = $currentTask.decoder
        TurboType = $currentTask.turboType
        CcsdsK = $currentTask.ccsdsK
        CcsdsRate = $currentTask.ccsdsRate
        StartSNR = $currentTask.startSNR
        EndSNR = $currentTask.endSNR
        Step = $currentTask.step
        Frames = $BatchSize
        Workers = 0
        OutputDir = $OutputDir
        SeedOffset = ($globalIndex * 10000) # Ensure each batch gets a vastly different seed base
    }
    
    # Execute parallel_sim
    $psArgs = @()
    foreach ($k in $simArgs.Keys) {
        $psArgs += "-$k"
        $psArgs += $simArgs[$k]
    }
    
    $p = Start-Process -FilePath "powershell" -ArgumentList "-ExecutionPolicy Bypass -File `"$ParallelSimPath`" $psArgs" -Wait -PassThru
    
    # Find the newly created merged CSV (parallel_sim creates one)
    $latestCsv = Get-ChildItem -Path $OutputDir -Filter "ber_*_merged.csv" | Sort-Object LastWriteTime -Descending | Select-Object -First 1
    
    if ($latestCsv -and $latestCsv.LastWriteTime -gt $batchStart) {
        Write-Host "  > Cycle done. Output: $($latestCsv.Name)"
        $state['LastOutputFile'] = $latestCsv.Name
    }
    
    $batchDuration = (Get-Date) - $batchStart
    $state['TotalDuration'] += $batchDuration
    $state['CompletedCycles']++
    $state['CompletedFrames'] += $thisBatch
    
    Write-Host "  > Task Cycle Completed." -ForegroundColor Green
}

# === 3. Generate Summary ===

$endTime = Get-Date
$totalDuration = $endTime - $startTime
$summaryFile = Join-Path $OutputDir ("overnight_summary_" + $startTime.ToString("yyyyMMdd_HHmmss") + ".md")

$summary = @()
$summary += "# Overnight Simulation Summary"
$summary += ""
$summary += "- **Start Time**: $($startTime)"
$summary += "- **End Time**: $($endTime)"
$summary += "- **Total Duration**: $($totalDuration.ToString())"
$summary += ""
$summary += "## Task Execution Report"
$summary += ""
$summary += "| Task ID | Name | Cycles | Status | Duration | Last Output |"
$summary += "|---|---|---|---|---|---|"

foreach ($t in $tasks) {
    $s = $taskState[$t.id]
    
    $row = "| $($t.id) | $($t.name) | $($s['CompletedCycles']) | $($s['Status']) | $($s['TotalDuration'].ToString()) | $($s['LastOutputFile']) |"
    $summary += $row
}

$summary | Set-Content $summaryFile

Write-Host "`n==========================================================" -ForegroundColor Cyan
Write-Host "   Overnight Batch Completed" -ForegroundColor Cyan
Write-Host "==========================================================" -ForegroundColor Cyan
Write-Host "Summary saved to: $summaryFile"

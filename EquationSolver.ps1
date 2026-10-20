param([switch]$Test)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

class QuantumEquationSolver {
    static [double] $Planck = 6.62607015e-34
    static [double] $HBar = 1.054571817e-34
    static [double] $ElectronVolt = 1.602176634e-19
    static [double] $SpeedOfLight = 299792458.0

    static [hashtable] SolveQuadratic([double]$a, [double]$b, [double]$c) {
        if ($a -eq 0) { throw 'Coefficient a must not be zero.' }
        $d = $b * $b - 4.0 * $a * $c
        if ($d -gt 0) {
            $s = [math]::Sqrt($d)
            return @{ Discriminant=$d; Roots=@((-$b+$s)/(2.0*$a),(-$b-$s)/(2.0*$a)); Kind='TwoRealRoots' }
        }
        if ($d -eq 0) {
            return @{ Discriminant=$d; Roots=@(-$b/(2.0*$a)); Kind='RepeatedRoot' }
        }
        $real = -$b/(2.0*$a)
        $imag = [math]::Sqrt(-$d)/(2.0*[math]::Abs($a))
        return @{ Discriminant=$d; Roots=@(@{Real=$real;Imaginary=$imag},@{Real=$real;Imaginary=(-$imag)}); Kind='ComplexConjugates' }
    }

    static [double] PhotonEnergyFromFrequency([double]$frequencyHz) {
        if ($frequencyHz -lt 0 -or [double]::IsNaN($frequencyHz) -or [double]::IsInfinity($frequencyHz)) { throw 'Frequency must be finite and non-negative.' }
        return [QuantumEquationSolver]::Planck * $frequencyHz
    }

    static [double] FrequencyFromPhotonEnergy([double]$energyJoules) {
        if ($energyJoules -lt 0 -or [double]::IsNaN($energyJoules) -or [double]::IsInfinity($energyJoules)) { throw 'Energy must be finite and non-negative.' }
        return $energyJoules / [QuantumEquationSolver]::Planck
    }

    static [double] DeBroglieWavelength([double]$momentumKgMPerS) {
        if ($momentumKgMPerS -le 0 -or [double]::IsNaN($momentumKgMPerS) -or [double]::IsInfinity($momentumKgMPerS)) { throw 'Momentum must be finite and positive.' }
        return [QuantumEquationSolver]::Planck / $momentumKgMPerS
    }

    static [double] PhotonEnergyFromWavelength([double]$wavelengthMeters) {
        if ($wavelengthMeters -le 0 -or [double]::IsNaN($wavelengthMeters) -or [double]::IsInfinity($wavelengthMeters)) { throw 'Wavelength must be finite and positive.' }
        return ([QuantumEquationSolver]::Planck * [QuantumEquationSolver]::SpeedOfLight) / $wavelengthMeters
    }

    static [double] ParticleInBoxEnergy([int]$n, [double]$massKg, [double]$lengthMeters) {
        if ($n -lt 1) { throw 'Quantum number n must be an integer >= 1.' }
        if ($massKg -le 0 -or [double]::IsNaN($massKg) -or [double]::IsInfinity($massKg)) { throw 'Mass must be finite and positive.' }
        if ($lengthMeters -le 0 -or [double]::IsNaN($lengthMeters) -or [double]::IsInfinity($lengthMeters)) { throw 'Box length must be finite and positive.' }
        return (($n * $n) * [QuantumEquationSolver]::Planck * [QuantumEquationSolver]::Planck) / (8.0 * $massKg * $lengthMeters * $lengthMeters)
    }

    static [double] HarmonicOscillatorEnergy([int]$n, [double]$angularFrequencyRadPerS) {
        if ($n -lt 0) { throw 'Quantum number n must be an integer >= 0.' }
        if ($angularFrequencyRadPerS -le 0 -or [double]::IsNaN($angularFrequencyRadPerS) -or [double]::IsInfinity($angularFrequencyRadPerS)) { throw 'Angular frequency must be finite and positive.' }
        return ($n + 0.5) * [QuantumEquationSolver]::HBar * $angularFrequencyRadPerS
    }
}

function Assert-Near([double]$Actual,[double]$Expected,[double]$Tolerance,[string]$Name) {
    if ([math]::Abs($Actual-$Expected) -gt $Tolerance) { throw "FAIL: $Name; actual=$Actual expected=$Expected" }
    Write-Host "PASS: $Name"
}
function Assert-Throws([scriptblock]$Action,[string]$Name) {
    $threw = $false
    try { & $Action | Out-Null } catch { $threw = $true }
    if (-not $threw) { throw "FAIL: $Name (expected exception)" }
    Write-Host "PASS: $Name"
}
function Invoke-Tests {
    $q = [QuantumEquationSolver]::SolveQuadratic(1,-5,6)
    Assert-Near ([double]$q.Roots[0]) 3 1e-12 'quadratic root 3'
    Assert-Near ([double]$q.Roots[1]) 2 1e-12 'quadratic root 2'
    $q = [QuantumEquationSolver]::SolveQuadratic(1,2,1)
    if ($q.Kind -ne 'RepeatedRoot') { throw 'FAIL: repeated root classification' }
    Assert-Near ([double]$q.Roots[0]) -1 1e-12 'repeated root -1'
    $q = [QuantumEquationSolver]::SolveQuadratic(1,0,1)
    if ($q.Kind -ne 'ComplexConjugates') { throw 'FAIL: complex root classification' }
    Assert-Near ([double]$q.Roots[0].Imaginary) 1 1e-12 'complex root +i'
    Assert-Near ([double]$q.Roots[1].Imaginary) -1 1e-12 'complex root -i'
    Assert-Near ([QuantumEquationSolver]::PhotonEnergyFromFrequency(1e14)) 6.62607015e-20 1e-32 'E=hf'
    Assert-Near ([QuantumEquationSolver]::FrequencyFromPhotonEnergy(6.62607015e-20)) 1e14 1 'f=E/h'
    Assert-Near ([QuantumEquationSolver]::DeBroglieWavelength(1e-24)) 6.62607015e-10 1e-22 'de Broglie wavelength'
    Assert-Near ([QuantumEquationSolver]::PhotonEnergyFromWavelength(QuantumEquationSolver::SpeedOfLight)) ([QuantumEquationSolver]::Planck) 1e-45 'E=hc/lambda at lambda=1m'
    Assert-Near ([QuantumEquationSolver]::ParticleInBoxEnergy(1,9.1093837e-31,1e-9)) 6.024e-20 1e-22 'particle in a box ground state (approx)'
    Assert-Near ([QuantumEquationSolver]::HarmonicOscillatorEnergy(0,1e14)) 5.272859085e-21 1e-33 'harmonic oscillator zero-point energy'
    Assert-Throws { [QuantumEquationSolver]::SolveQuadratic(0,2,1) } 'reject invalid quadratic'
    Assert-Throws { [QuantumEquationSolver]::DeBroglieWavelength(0) } 'reject zero momentum'
    Assert-Throws { [QuantumEquationSolver]::ParticleInBoxEnergy(0,1,1) } 'reject invalid quantum number'
    Write-Host 'ALL TESTS PASSED' -ForegroundColor Green
}

if ($Test) {
    Invoke-Tests
    exit 0
}

Write-Host ''
Write-Host 'HEAVENZFIRE | QUANTUM EQUATION SOLVER' -ForegroundColor Cyan
Write-Host 'Numerical educational tool — SI units unless stated otherwise.'
Write-Host '1) Solve a*x^2 + b*x + c = 0'
Write-Host '2) Photon energy from frequency (E = h*f)'
Write-Host '3) Photon energy from wavelength (E = h*c/lambda)'
Write-Host '4) de Broglie wavelength (lambda = h/p)'
Write-Host '5) Particle in a 1D infinite potential well'
Write-Host '6) Quantum harmonic oscillator energy level'
$choice = Read-Host 'Select 1-6'
$culture = [Globalization.CultureInfo]::InvariantCulture
switch ($choice) {
    '1' {
        $a = [double]::Parse((Read-Host 'a'),$culture); $b = [double]::Parse((Read-Host 'b'),$culture); $c = [double]::Parse((Read-Host 'c'),$culture)
        $r = [QuantumEquationSolver]::SolveQuadratic($a,$b,$c)
        Write-Host ("Type: {0}; discriminant: {1}" -f $r.Kind,$r.Discriminant)
        $r.Roots | ForEach-Object { if ($_ -is [hashtable]) { Write-Host ("x = {0} {1} {2}i" -f $_.Real,$(if($_.Imaginary -ge 0){'+'}else{''}),$_.Imaginary) } else { Write-Host ("x = {0}" -f $_) } }
    }
    '2' { $f = [double]::Parse((Read-Host 'Frequency in Hz'),$culture); $e=[QuantumEquationSolver]::PhotonEnergyFromFrequency($f); Write-Host ("Energy: {0} J ({1} eV)" -f $e,($e/[QuantumEquationSolver]::ElectronVolt)) }
    '3' { $l = [double]::Parse((Read-Host 'Wavelength in meters'),$culture); $e=[QuantumEquationSolver]::PhotonEnergyFromWavelength($l); Write-Host ("Energy: {0} J ({1} eV)" -f $e,($e/[QuantumEquationSolver]::ElectronVolt)) }
    '4' { $p = [double]::Parse((Read-Host 'Momentum in kg*m/s'),$culture); Write-Host ("de Broglie wavelength: {0} m" -f [QuantumEquationSolver]::DeBroglieWavelength($p)) }
    '5' { $n = [int]::Parse((Read-Host 'Quantum number n (integer >= 1)'),$culture); $m=[double]::Parse((Read-Host 'Particle mass in kg'),$culture); $l=[double]::Parse((Read-Host 'Box length in meters'),$culture); $e=[QuantumEquationSolver]::ParticleInBoxEnergy($n,$m,$l); Write-Host ("E_n = {0} J ({1} eV)" -f $e,($e/[QuantumEquationSolver]::ElectronVolt)) }
    '6' { $n = [int]::Parse((Read-Host 'Quantum number n (integer >= 0)'),$culture); $w=[double]::Parse((Read-Host 'Angular frequency in rad/s'),$culture); $e=[QuantumEquationSolver]::HarmonicOscillatorEnergy($n,$w); Write-Host ("E_n = {0} J ({1} eV)" -f $e,($e/[QuantumEquationSolver]::ElectronVolt)) }
    default { throw 'Choose a menu option from 1 to 6.' }
}

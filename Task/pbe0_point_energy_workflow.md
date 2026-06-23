# PBE0 Point Energy Calculation Workflow

## 1. Prepare Input Structures

Put each input structure in its own folder, ideally as `CONTCAR`.

```text
source_root/
  SYSTEM_A/CONTCAR
  SYSTEM_B/CONTCAR
```

If the input is a `POSCAR`, copy or rename it to `CONTCAR` for the current generator script.

## 2. Generate VASP Inputs With CSPToolbox

Use:

```bash
/opt/anaconda3/bin/python /Users/zianzhan/Desktop/CSP_sandbox/CSPToolbox/Source/CLI_scripts/generate_pbe0_vasp_inputs_from_latest_contcar.py \
  --source-root /path/to/source_root \
  --output-root /path/to/02_PBE0 \
  --incar-template /Users/zianzhan/Desktop/CSP_sandbox/CSPToolbox/Template/VASP_input/PBE0_INCAR_1000eV_point_calc \
  --kpoint-density 0.05
```

`--kpoint-density 0.05` means a k-point spacing of `0.05 * 2pi A^-1`.

## 3. Check Generated Files

Each generated system folder should contain:

```text
POSCAR
INCAR
KPOINTS
POTCAR
run_vasp_1000eV.sh
```

CSPToolbox writes `POSCAR` in Cartesian coordinates. This is fine; Direct vs Cartesian should not change the calculation if the positions are equivalent.

## 4. Main PBE0 INCAR Settings

The point-energy setup uses:

```text
ENCUT = 1000
NSW = 0
IBRION = -1
LHFCALC = .TRUE.
AEXX = 0.25
HFSCREEN = 0.0
IVDW = 12
LWAVE = .TRUE.
LCHARG = .TRUE.
```

For IBZKPT symmetry errors, add:

```text
ISYM = 0
```

## 5. Check Parallel Settings

Typical setup used here:

```text
#$ -pe mpi 80
#$ -l mem=4G
```

For harder reruns, `mpi 160` was used. For insufficient-memory cases, set:

```text
KPAR = 1
```

## 6. Add STOPCAR Walltime Safeguard

For a `24:00:00` walltime job, add logic to write `STOPCAR` before walltime is hit:

```text
LABORT = .TRUE.
```

The current convention is to trigger this after `sleep 22h`.

## 7. Upload To Young

Typical pattern:

```bash
tar -czf 02_PBE0.tar.gz 02_PBE0
scp 02_PBE0.tar.gz y:/home/mmm1598/Scratch/
ssh y
cd /home/mmm1598/Scratch
tar -xzf 02_PBE0.tar.gz
```

Be careful not to overwrite existing completed result folders.

## 8. Submit Jobs

From the remote calculation root:

```bash
cd /home/mmm1598/Scratch/02_PBE0
qsub SYSTEM/run_vasp_1000eV.sh
```

For multiple jobs, use a batch script targeting only the intended systems.

## 9. Check Completion

A clean completed calculation should have:

- `OUTCAR` ending with `General timing and accounting`
- final `F=` energy in `vasp.out`
- no fatal marker such as `IBZKPT`, insufficient memory, killed, etc.

The final PBE0 energy is taken from the last `F=` line in `vasp.out`.

## 10. Postprocess Energy

For per-molecule or per-`Z` lattice energy:

```text
E_latt_PBE0_kJ/mol = (PBE0_eV / Z) * 96.485
```

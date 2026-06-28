# SCR ensemble run on HPC

This pipeline is independent of `run-producer_historic.py`, `run-consumer.py`,
and `sim_setups_historic.csv`. It uses the dedicated `sim_setups_scr.csv`,
`sim_scr.json`, and `crop_scr.json` files.
It creates 50 representative SCR/BKR points, applies one LHS parameter group to
each point, runs the 1980-2020 climate period, and writes 40 complete winter-wheat
seasons (growing years 1981-2020) per scenario CSV.

For each polygon, the representative-point builder first keeps cells with valid
BUEK soil and crop-mask value 1. It finds the polygon's modal BUEK soil id and
selects the cell with that soil id nearest to the valid-cell centroid. The nearest
DWD climate cell is assigned later by the same mapping used by the historic
producer.

## 1. Create deterministic inputs

Run this once in the project directory:

```bash
python3 prepare_scr_inputs.py parameter-groups=20 lhs-seed=42
```

Use `parameter-groups=25` or `parameter-groups=30` for a larger ensemble. Use a
new output directory whenever the parameter-group count, seed, templates, or
simulation dates change.

The default 20-group design is balanced between 10 rainfed and 10 irrigated
groups. Irrigated groups apply 4-8 mm per automatic-irrigation event with a
soil-moisture threshold of 0.075-0.14. Daily `Irrig` output records the amount
actually applied by MONICA.

## 2. Preflight on the HPC filesystem

The full preflight checks all representative soil profiles and every selected
climate file, including its date coverage. It does not send MONICA jobs.

```bash
python3 run-producer_scr.py \
  mode=remoteProducer-remoteMonica \
  out=scr-csv-out-20 \
  check-climate-files=true \
  preflight-only=true
```

Do not start the run if this command reports a missing climate file, missing soil
profile, or date-coverage error.

## 3. End-to-end smoke test

Before the full ensemble, run one complete 1980-2020 scenario through the real
HPC proxies and MONICA worker. Start the consumer first:

```bash
python3 run-consumer_scr.py \
  server=login01.cluster.zalf.de port=7777 shared_id=scr-lhs-smoke \
  out=scr-smoke max-points=1 max-params=1 timeout=7200000
```

Then send the matching single environment:

```bash
python3 run-producer_scr.py \
  mode=remoteProducer-remoteMonica \
  server=login01.cluster.zalf.de server-port=6666 shared_id=scr-lhs-smoke \
  out=scr-smoke max-points=1 max-params=1 check-climate-files=true
```

Continue only when the consumer exits with code 0 and
`scr-smoke/BKR_101/BKR101_P001_daily.csv` contains 40 growing years.

## 4. Full run

Start the MONICA proxies/workers as usual. Start exactly one consumer first:

```bash
python3 run-consumer_scr.py \
  server=login01.cluster.zalf.de \
  port=7777 \
  shared_id=scr-lhs \
  out=scr-csv-out-20 \
  timeout=7200000
```

Then start the producer in another HPC shell/job:

```bash
python3 run-producer_scr.py \
  mode=remoteProducer-remoteMonica \
  server=login01.cluster.zalf.de \
  server-port=6666 \
  shared_id=scr-lhs \
  out=scr-csv-out-20 \
  check-climate-files=true
```

The producer and consumer must use the same `points-file`, `params-file`,
`max-points`, `max-params`, `shared_id`, `out`, and `expected-seasons` values.

## 5. Resume after interruption

Run the same consumer and producer commands again. The producer skips final CSV
files already present. The consumer writes each scenario to `.tmp`, validates 40
complete seasons, calls `fsync`, and only then atomically renames it to the final
CSV. A stale `.tmp` file is therefore never treated as complete.

The output layout is:

```text
scr-csv-out-20/
  run_manifest.json
  run_status.csv
  missing_scenarios.csv
  BKR_101/BKR101_P001_daily.csv
  BKR_101/BKR101_P002_daily.csv
  ...
```

With 20 parameter groups, the run contains 1,000 MONICA environments and 40,000
complete crop-season samples. A non-empty `missing_scenarios.csv` or a non-zero
consumer exit code means the run is not complete.

## Local flat climate directory

MONICA reads `.csv.gz` directly; decompression is not required. For climate files
stored flat under `data/dwd/dwd`, use:

```text
mode=localProducer-localMonica flat-climate-dir=data/dwd/dwd server=localhost
```

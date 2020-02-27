cat input/dt_* > input/dt.cc
cat input/event_* > input/event.dat
hypoDD hypoDD.inp
python reloc2csv.py
mv hypoDD.log output

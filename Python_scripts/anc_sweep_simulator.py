#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Usage: anc_sweep_simulator.py -n <INT> -m <FLT> -r <FLT> -l <INT> -s <STR> -i <INT> [-h]

  [Options]
    -n, --Ne <INT>                             Effective population size
    -m, --mu <FLT>                             Mutation rate
    -r, --recomb_rate <FLT>                    Recombination rate
    -l, --seq_length <INT>                     Sequence length
    -s, --sweep_params <STR>                   Sweep parameters: s,Ta
    -i, --replicates <INT>                     The number of times to repeat the simulation
    -h, --help                                 Show this message

"""

from docopt import docopt
import msprime
import sys

if __name__ == '__main__':
  __version__ = '0.1'
  args = docopt(__doc__)
  #print(args)
  Ne = int(args["--Ne"])
  mu = float(args["--mu"])
  seq_length = int(args["--seq_length"])
  recomb_rate = float(args["--recomb_rate"])
  sel, Ta = [float(k) for k in args["--sweep_params"].split(",")]
  replicates = int(args["--replicates"])
  demography = msprime.Demography()
  demography.add_population(name="X", initial_size=Ne)
  #print(demography.debug())
  if sel == 0:
    ts_reps = msprime.sim_ancestry(samples={"X": 4}, demography=demography, 
    sequence_length=seq_length, recombination_rate=recomb_rate, 
    model=[msprime.StandardCoalescent()], num_replicates=replicates)
  else:
    sweep_model = msprime.SweepGenicSelection(position = seq_length / 2, 
    start_frequency = 1.0 / (2 * Ne), 
    end_frequency = 1.0 - (1.0 / (2 * Ne)), 
    s=sel, 
    dt=1e-6)
    ts_reps = msprime.sim_ancestry(samples={"X": 4}, demography=demography, 
    sequence_length=seq_length, recombination_rate=recomb_rate, 
    model=[msprime.StandardCoalescent(duration=Ta), 
    sweep_model, 
    msprime.StandardCoalescent()], num_replicates=replicates)
  for i, ts in enumerate(ts_reps):
    mutated_ts = msprime.sim_mutations(ts, rate=mu)
    mutated_ts.write_vcf(sys.stdout, contig_id=str(i))

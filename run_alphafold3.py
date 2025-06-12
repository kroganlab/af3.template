#!/usr/bin/python3
#
#$ -S /usr/bin/python3
#$ -q gpu.q
#$ -N alphafold 
#$ -cwd
#$ -l h_rt=2:00:00
#$ -l mem_free=64G
#$ -l scratch=50G
#$ -l compute_cap=80,gpu_mem=40G
###$ -pe smp 2
#
# Compute cap for A100 GPU is 8.0 (40 or 80 GB), for A40 GPU is 8.6 (48 GB).
#
# The "-pe smp 4" means allocate 4 slots.  On a CPU queue that means 4 cores,
# but on the gpu.q in means 4 GPUs.
#
# mem_free is per-slot.
#
# Adapted from alphafold/docker/run_alphafold.py script.
# Original version runs AlphaFold using a docker image.
# This adapted version uses a singularity image with defaults
# set for the UCSF Wynton cluster.
#

"""Singularity launch script for Alphafold."""

def parse_args():
  import argparse

  parser = argparse.ArgumentParser(description='Run AlphaFold structure prediction using singularity image.')
  ###
  # added params
  ###
  parser.add_argument(
    '--job_id', type = int, required = True,
    help = 'job id to look up sequence names in AlphaFoldJobList')

  parser.add_argument(
    '--master_fasta', required=True, default = "masterFasta.fasta",
    help = "fasta file with sequence info relevant to names in AlphaFoldJobList")
    
  parser.add_argument(
    '--jobTable', required=True, default= "AlphaFoldJobList.csv",
    help = "path to csv formatte file with columns ID, seq1.name, seq2.name, seq3.name etc..")  
  
  parser.add_argument(
    '--alignmentRepo', default  =  "/wynton/group/krogan/mgordon/data/AF3_A3M_MSAs",
    help = "path to directory containing previously generated MSAs in A3M format")
  
  parser.add_argument(
    '--dependenciesRepo', default  =  "/wynton/group/krogan/mgordon/data",
    help = "path repo containing dependencies for AF3 pipeline")
  
  parser.add_argument(
    '--nSeeds', default  =  5,
    help = "Number of random seeds for initalizations for AF model")
  
  parser.add_argument(
    '--json_template', default  =  "/wynton/group/krogan/mgordon/data/af3_template.json",
    help = "JSON template to be populated for AF runs")  

  parser.add_argument(
    '--blockChainMSA', type = bool, required = False, default = False,
    help = 'Add padding to the unpaired MSAs (By default AF3 concatenates unpaired MSAs from different chains).\nNote: this will reduce unpaired MSA depth by 1/N chains')
    
  # dont think this is needed anymore...  
  # maybe reuse, should just look for the input json in the outDir
  parser.add_argument(
    '--setup_job', type = str_to_bool, default = True,
    help = "Disable to skip the setup steps and run inference only")
  
  ###
  # end added params
  ###

  # Must specify either json_path or input_dir
  parser.add_argument(
    '--json_path',
    help='Paths to the input JSON file')

  parser.add_argument(
    '--input_dir',
    help='Paths to the directory containing input JSON files')

  parser.add_argument(
    '--output_dir',
    default = '.',
    help='Paths to a directory where the results will be saved')

  from os.path import expanduser, isdir
  # AF3 weights taken from /wynton/home/ferrin/goddard/af3_weights
  parser.add_argument(
    '--model_dir',
    default = '/wynton/home/ferrin/goddard/af3_weights',
    help='Path to the model to use for inference.')

  parser.add_argument(
    '--flash_attention_implementation',
    default='triton',
    choices=['triton', 'cudnn', 'xla'],
    help=(
        "Flash attention implementation to use. 'triton' and 'cudnn' uses a"
        ' Triton and cuDNN flash attention implementation, respectively. The'
        ' Triton kernel is fastest and has been tested more thoroughly. The'
        " Triton and cuDNN kernels require Ampere GPUs or later. 'xla' uses an"
        ' XLA attention implementation (no flash attention) and is portable'
        ' across GPU devices.'
    ),
  )

  # Control which stages to run.
  parser.add_argument(
    '--run_data_pipeline',
    default=True, type=str_to_bool, 
    help='Whether to run the data pipeline on the fold inputs.',
  )

  parser.add_argument(
    '--run_inference',
    default=True, type=str_to_bool, 
    help='Whether to run inference on the fold inputs.',
  )

  parser.add_argument(
    '--db_dir',
    default = '/wynton/group/databases/alphafold3',
    help = 'Path to the directory containing the databases.',
  )

  # Number of CPUs to use for MSA tools.
  import multiprocessing
  parser.add_argument(
    '--jackhmmer_n_cpu',
    default = min(multiprocessing.cpu_count(), 8),
    help = 'Number of CPUs to use for Jackhmmer. Default to min(cpu_count, 8). Going'
    ' beyond 8 CPUs provides very little additional speedup.',
  )

  parser.add_argument(
    '--nhmmer_n_cpu',
    default = min(multiprocessing.cpu_count(), 8),
    help = 'Number of CPUs to use for Nhmmer. Default to min(cpu_count, 8). Going'
           ' beyond 8 CPUs provides very little additional speedup.',
  )

  # Compilation cache
  parser.add_argument(
    '--jax_compilation_cache_dir',
    default = None,
    help ='Path to a directory for the JAX compilation cache.',
  )

  import os
  parser.add_argument(
    '--gpu_devices', default=os.environ.get('SGE_GPU', '0'),
    help='Comma separated list GPU identifiers to set environment variable CUDA_VISIBLE_DEVICES.')

  parser.add_argument(
    '--singularity_image_path',
    help='Path to the AlphaFold singularity image.')

  parser.add_argument(
    '--use_a100_80gb_settings', type=str_to_bool,
    help='Use AlphaFold 3 settings for A100 80 GB graphics.  If not set use A100 40 GB settings.')
  
  args = parser.parse_args()
  return args

def str_to_bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        import argparse
        raise argparse.ArgumentTypeError('Boolean value expected.')

def main():
  import os
  args = parse_args()

  if args.model_dir is None:
    raise RuntimeError("If you do not have the AlphaFold 3 neural net weights (also called model parameters) in directory ~/af3_weights then you must provide the --model_dir option to specify the location of the weights.  To obtain the weights you must submit a request to Google https://forms.gle/svvpY4u2jsHEwWYS6 as described on the AlphaFold 3 Github page https://github.com/google-deepmind/alphafold3.")

  # MG functions to stage AF runs
  import Alphafold3_utils

  args.json_path = Alphafold3_utils.af3_setupJob(job_id=args.job_id, 
                                                 setup_job=args.setup_job,
                                                 jobTable=args.jobTable,
                                                 master_fasta=args.master_fasta,
                                                 alignmentRepo=args.alignmentRepo,
                                                 nSeeds=args.nSeeds,
                                                 json_template=args.json_template,
                                                 output_dir=args.output_dir,
                                                 blockChainMSA=args.blockChainMSA
                                                 )
  
  print(f'Running AF on following input:\n{args.json_path}')
  
  mounts = []
  command_args = []

  run_args = ['json_path', 'input_dir', 'db_dir', 'output_dir', 'model_dir', 'flash_attention_implementation',
              'run_data_pipeline', 'run_inference', 'jackhmmer_n_cpu', 'nhmmer_n_cpu', 'jax_compilation_cache_dir']
  for arg_name in run_args:
      if getattr(args, arg_name) is not None:
          command_args.append(f'--{arg_name}={getattr(args,arg_name)}')

  if args.json_path is None and args.input_dir is None:
     raise ValueError('Exactly one of --json_path or --input_dir must be specified.')

  env_vars = {
          'CUDA_VISIBLE_DEVICES': args.gpu_devices,
          'NVIDIA_VISIBLE_DEVICES': args.gpu_devices,
  }
  env_vals = ','.join('%s=%s' % (key,value) for key,value in env_vars.items())

  # AlphaFold uses Python tempfile which uses TMPDIR env variable
  # which is /scratch/job-id-string on wynton.  Otherwise Python will use /tmp
  # which is only 4-8 GB on wynton and will cause write errors on large sequences.
  import os
  tempdir = os.environ.get('TMPDIR', '/scratch')

  # Mount AlphaFold databases, models directory, current directory, scratch directory
  bind_directories = [args.db_dir, args.model_dir, args.alignmentRepo, args.dependenciesRepo, os.getcwd(), tempdir]

  # Bind parent directories of fasta file locations, and output directory..
  from os.path import isabs, dirname
  if args.json_path and isabs(args.json_path) and dirname(args.json_path):
      bind_directories.append(dirname(args.json_path))
  if args.input_dir and isabs(args.input_dir):
    bind_directories.append(args.input_dir)
  if args.output_dir and isabs(args.output_dir):
    bind_directories.append(dirname(args.output_dir))
  if args.jax_compilation_cache_dir and isabs(args.jax_compilation_cache_dir):
    bind_directories.append(args.jax_compilation_cache_dir)

  if args.singularity_image_path:
    singularity_image_path = args.singularity_image_path
  elif args.use_a100_80gb_settings:
    singularity_image_path = '/wynton/home/ferrin/goddard/alphafold_singularity/alphafold3_80gb.sif'
  else:
    singularity_image_path = '/wynton/group/krogan/mgordon/projects/112624_MGordon_AF3_pipeline/testing/af3_env/af3_env_mg.sif'
    #singularity_image_path = '/wynton/group/krogan/mgordon/projects/112624_MGordon_AF3_pipeline/testing/af3_env/af3_env_mg_dev'
    #singularity_image_path = '/wynton/group/krogan/mgordon/projects/112624_MGordon_AF3_pipeline/script/af3_40gb_dec_4_2024_sandbox'

  subprocess_args = ['singularity',
          'exec',
          '--nv',  # Use Nvidia container library to use CUDA
          '-B "%s"' % ','.join(bind_directories),
          '--env %s' % env_vals,
          singularity_image_path,
          'python', '/app/alphafold/run_alphafold.py',
        ] + command_args
  cmd = ' '.join(subprocess_args)

  print(f"\nAF3 command executed:\n{cmd}\n")

  from subprocess import run
  import sys
  run('module load cuda/12.2 ; %s' % cmd,
      stdout = sys.stdout, stderr = sys.stderr,
      shell = True,  # module command is a csh alias on Wynton
      executable = '/bin/csh',
      check = True)

  seqlist = os.path.basename(args.json_path).split(".af3_input.json")[0].split("__")
  runOutdir = Alphafold3_utils.af3_alphaFoldRunOutputDirectory(seqlist, args.output_dir)

  # for now only enable MSA recovery if blockChain unset (dont want to include these padded MSA in our repo....)
  # TODO move the recovery function to before AF inference,
  #  may need to modify source code...
  if args.blockChainMSA:
    print(f"blockChainMSA is enabled. Not extracting MSAs..")
  else:
    print(f"Recovering MSAs and saving to {args.alignmentRepo}...")
    Alphafold3_utils.af3_captureMSAs(output_dir=runOutdir, 
                                     alignmentRepo=args.alignmentRepo)
  
  print('')
  print('Extracting model PTM & iPTM scores...')
  Alphafold3_utils.af3_captureSummaryScores(output_dir=runOutdir)
  print('')
  print('Generating MSA & PAE plots for top ranking model...')
  Alphafold3_utils.generate_MSAandPAEplots(outDir=runOutdir)
  print('')
  print('Getting inter-chain contacts and PAE...')
  Alphafold3_utils.get_interchainContactsPAE(outDir=runOutdir)
  print('')
  print('plotting inter-chain distances...')
  Alphafold3_utils.plot_interChainDistances(outDir=runOutdir)

if __name__ == '__main__':
  main()
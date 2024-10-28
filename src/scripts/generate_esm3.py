import sys
import os
os.environ["TOKENIZERS_PARALLELISM"] = "false"
# os.environ["TRANSFORMERS_CACHE"] = "/nethome/kgeorge/workspace/DomainPrediction/Data/esm3_experiments/mi_exp"
# os.environ["MPLCONFIGDIR"] = "/nethome/kgeorge/workspace/DomainPrediction/Data/esm3_experiments/mi_exp"
sys.path.append('../../esm') ## ignore if you are installing esm3 and use huggingface_hub login()
sys.path.append('..')

import numpy as np
import torch

from esm.models.esm3 import ESM3
from esm.sdk.api import ESM3InferenceClient, ESMProtein, GenerationConfig
from esm.utils.structure.protein_chain import ProteinChain

from utils import helper

protein = ProteinChain.from_pdb('../../Data/1bag.protein.pdb')

sequence_prompt = ''.join(['_'  for _ in range(len(protein))]) ## fully masked prompt
structure_prompt = torch.tensor(protein.atom37_positions) ## get structure as input

model: ESM3InferenceClient = ESM3.from_pretrained("esm3_sm_open_v1").to("cuda")

fasta_file = '../../Data/esm3_gen.fasta' ## file loc

N_GENERATIONS = 200
temperature = 0.5
for idx in range(N_GENERATIONS):
    sequence_prediction_config = GenerationConfig(
        track="sequence", 
        num_steps=sequence_prompt.count("_") // 1, 
        temperature=temperature
    )
    esm_protein = ESMProtein(sequence=sequence_prompt, coordinates=structure_prompt)
    generated_protein = model.generate(esm_protein, sequence_prediction_config)

    seq_dict = {}
    gen_idx = f'amylase_esm3_temp_{temperature}_gen_{idx}'
    seq_dict[gen_idx] = generated_protein.sequence

    print(gen_idx)

    helper.create_fasta(seq_dict, fasta_file, append=True)


import sys
sys.path.append('..')

import os
from tqdm import tqdm
import numpy as np
import torch
import esm

from utils import helper

class ESM2():
    def __init__(self, model_path, device='cpu') -> None:
        # self.model, self.alphabet = esm.pretrained.load_model_and_alphabet(model_path)
        self.model, self.alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        self.batch_converter = self.alphabet.get_batch_converter()
        self.model.eval()
        self.device = device

        if self.device == 'gpu':
            self.model.cuda()

        self.tok_to_idx = self.alphabet.tok_to_idx
        self.idx_to_tok = {v:k for k,v in self.tok_to_idx.items()}

    def get_res(self, sequence):
        data = [
            ("protein1", sequence)
        ]
        batch_labels, batch_strs, batch_tokens = self.batch_converter(data)
        batch_lens = (batch_tokens != self.alphabet.padding_idx).sum(1)

        if self.device == 'gpu':
            batch_tokens = batch_tokens.cuda()

        with torch.no_grad():
            results = self.model(batch_tokens, repr_layers=[33], return_contacts=True)

        return results

    def get_logits(self, sequence):

        results = self.get_res(sequence)
        return results['logits']

    def get_prob(self, sequence):
        logits = self.get_logits(sequence)
        prob = torch.nn.functional.softmax(logits, dim=-1)[0, 1:-1, :] # 1st and last are start and end tokens

        return prob.cpu().numpy()
    
def compute_perplexity(model, sequence, mask_token='<mask>'):
    '''
        pseudoperplexity(x) = exp( -1/L \sum_{i=1}_{L} [log( p(x_{i}|x_{j!=i}) )] )
    '''
    
    sum_log = 0
    for pos in tqdm(range(len(sequence))):
        masked_query = list(sequence)
        assert mask_token not in masked_query
        masked_query[pos] = mask_token
        masked_query = ''.join(masked_query)
        prob = model.get_prob(sequence=masked_query)

        assert prob.shape[0] == len(sequence)

        prob_pos = np.log(prob[pos, model.tok_to_idx[sequence[pos]]])
        
        sum_log += prob_pos

    return np.exp(-1*sum_log/len(sequence))

model_path = '/data/users/kgeorge/workspace/esm2/checkpoints/esm2_t33_650M_UR50D.pt'
esm2 = ESM2(model_path = model_path, device='gpu') # model_path dosen't matter - I have locally stored esm2

fasta_path = '../../Data/sequences_425.fasta'
meta_file = '../../Data/sequences_425.metadata.json'

records = helper.read_fasta(fasta_path)
for rec in records:
    perplexity = compute_perplexity(esm2, str(rec.seq))

    print(rec.id, perplexity)
    helper.update_metadata_json(meta_file, rec.id, 'esm2_650M_perplexity', perplexity, force=False)
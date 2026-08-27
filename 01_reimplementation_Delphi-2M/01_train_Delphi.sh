#!/bin/bash

cd "<path_to_file>"

python train.py config/train_delphi.py --device=cuda --out_dir=Delphi-2M-real
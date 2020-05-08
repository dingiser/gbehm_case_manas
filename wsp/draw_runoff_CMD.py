
from GBEHM.draw import draw_runoff
import argparse

parser=argparse.ArgumentParser("draw_runoff_CMD.py")
parser.add_argument("simu",help="simulated runoff by GBEHM")
parser.add_argument('-v',"--val",required=False,help="validation runoff")
parser.add_argument("-s",'--save',help='output file name.')
args=parser.parse_args()

draw_runoff(args.simu,val_runoff=args.val,save=args.save)
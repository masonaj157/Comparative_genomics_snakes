#!/usr/bin/env python

import argparse
import csv


parser = argparse.ArgumentParser(description='Builds the ctl file for PAML to help automate CodeML analyses')
parser.add_argument("-s","--seqfile",
					type=str,
					default='',
					help="Sequence file in phylip format")
parser.add_argument("-t","--treefile",
					type=str,
					default='',
					help="Tree file for the genetree. Newick Format")
parser.add_argument("-o","--outfile",
					type=str,
					default='output',
					help="name of the main result file (e.x. output)")
parser.add_argument("-n","--noisy",
					type=str,
					default='9',
					help="how much stuff you want on the screen")
parser.add_argument("-v","--verbose",
					type=str,
					default='1',
					help="from codeml ctl - 0: concise; 1: detailed, 2: too much")
parser.add_argument("-r","--runmode",
					type=str,
					default='0',
					help="0: user tree;  1: semi-automatic;  2: automatic;  3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise")
parser.add_argument("-st","--seqtype",
					type=str,
					default='1',
					help="1:codons; 2:AAs; 3:codons-->AAs")
parser.add_argument("-cf","--codonfreq",
					type=str,
					default='2',
					help="0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table")
parser.add_argument("-nd","--ndata",
					type=str,
					default='1',
					help="number of datasets to analyze")
parser.add_argument("-c","--clock",
					type=str,
					default='0',
					help="0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis")
parser.add_argument("-ad","--aaDist",
					type=str,
					default='0',
					help="0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a")
parser.add_argument("-ar","--aaRatefile",
					type=str,
					default='0',
					help="0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a")
parser.add_argument("-m","--model",
					type=str,
					default='2',
					help=" models for codons: 0:one, 1:b, 2:2 or more dN/dS ratios for branches | models for AAs or codon-translated AAs: 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F, 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)")
parser.add_argument("-ns","--NSsites",
					type=str,
					default='0',
					help="0:one w;1:neutral;2:selection; 3:discrete;4:freqs; 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma; 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1; 13:3normal>0")
parser.add_argument("-i","--icode",
					type=str,
					default='0',
					help=" Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt., 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 10: blepharisma nu., These codes correspond to transl_table 1 to 11 of GENEBANK.")
parser.add_argument("-mg","--Mgene",
					type=str,
					default='0',
					help=" codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff | AA: 0:rates, 1:separate")
parser.add_argument("-fk","--fix_kappa",
					type=str,
					default='0',
					help="1: kappa fixed, 0: kappa to be estimated")
parser.add_argument("-k","--kappa",
					type=str,
					default='2',
					help="initial or fixed kappa")
parser.add_argument("-fo","--fix_omega",
					type=str,
					default='0',
					help="1: omega or omega_1 fixed, 0: estimate ")
parser.add_argument("-om","--omega",
					type=str,
					default='0.4',
					help="initial or fixed omega, for codons or codon-based AAs")
parser.add_argument("-fa","--fix_alpha",
					type=str,
					default='1',
					help="0: estimate gamma shape parameter; 1: fix it at alpha")
parser.add_argument("-a","--alpha",
					type=str,
					default='0',
					help="initial or fixed alpha, 0:infinity (constant rate)")
parser.add_argument("-ma","--Malpha",
					type=str,
					default='0',
					help="different alphas for genes")
parser.add_argument("-ncg","--ncatG",
					type=str,
					default='8',
					help="# of categories in dG of NSsites models")
parser.add_argument("-gse","--getSE",
					type=str,
					default='0',
					help="0: don't want them, 1: want S.E.s of estimates")
parser.add_argument("-ra","--rateAncestor",
					type=str,
					default='1',
					help="(0,1,2): rates (alpha>0) or ancestral states (1 or 2)")
parser.add_argument("-sd","--Small_Diff",
					type=str,
					default='0.5e-6',
					help="...")
parser.add_argument("-cd","--cleandata",
					type=str,
					default='1',
					help="remove sites with ambiguity data (1:yes, 0:no)")
parser.add_argument("-met","--method",
					type=str,
					default='0',
					help="Optimization method 0: simultaneous; 1: one branch a time")																				
args=parser.parse_args()

## Check arguments for errors
if args.seqfile == '':
	raise ValueError('No seqfile specified')
	quit()


if args.treefile == '':
	raise ValueError('No treefile specified')
	quit()


if args.noisy not in ['0','1','2','3','9']:
	raise ValueError('unrecognized value for noisy')
	quit()
	

if args.verbose not in ['0','1','2']:
	raise ValueError('unrecognized value for verbose')
	quit()


if args.runmode not in ['0','1','2','3','4','5','-2']:
	raise ValueError('unrecognized value for runmode')
	quit()


if args.seqtype not in ['1','2','3']:
	raise ValueError('unrecognized value for seqtype')
	quit()
	

if args.codonfreq not in ['0','1','2','3']:
	raise ValueError('unrecognized value for CodonFreq')
	quit()

if args.clock not in ['0','1','2','3']:
	raise ValueError('unrecognized value for clock')
	quit()

if args.aaDist not in ['0','+','-','1','2','3','4','5','6']:
	raise ValueError('unrecognized value for clock')
	quit()


if args.model not in ['0','1','2','3','6','7','8','9']:
	raise ValueError('unrecognized value for model')
	quit()

if args.NSsites not in ['0','1','2','3','4','5','6','7','8','9','10','11','12','13']:
	raise ValueError('unrecognized value for NSites')
	quit()

if args.icode not in ['0','1','2','3','4','5','6','7','8','9','10']:
	raise ValueError('unrecognized value for icode')
	quit()


if args.Mgene not in ['0','1','2','3','4']:
	raise ValueError('unrecognized value for Mgene')
	quit()


if args.fix_kappa not in ['0','1']:
	raise ValueError('unrecognized value for fix_kappa')
	quit()
	

if args.fix_omega not in ['0','1']:
	raise ValueError('unrecognized value for fix_omega')
	quit()


if args.fix_omega not in ['0','1']:
	raise ValueError('unrecognized value for fix_omega')
	quit()


if args.fix_alpha not in ['0','1']:
	raise ValueError('unrecognized value for fix_alpha')
	quit()


if args.getSE not in ['0','1']:
	raise ValueError('unrecognized value for getSE')
	quit()


if args.rateAncestor not in ['0','1','2']:
	raise ValueError('unrecognized value for rateAncestor')
	quit()


if args.cleandata not in ['0','1']:
	raise ValueError('unrecognized value for cleandata')
	quit()


if args.method not in ['0','1']:
	raise ValueError('unrecognized value for cleandata')
	quit()
	


## Now write the codemlfile
outFile = open('codeml.ctl', "w")
outFile.write('      seqfile = ' + args.seqfile + '\n')
outFile.write('     treefile = ' + args.treefile + '\n')
outFile.write('      outfile = ' + args.outfile + '\n')
outFile.write('\n')
outFile.write('        noisy = ' + args.noisy + '\n')
outFile.write('      verbose = ' + args.verbose + '\n')
outFile.write('      runmode = ' + args.runmode + '\n')
outFile.write('\n')
outFile.write('\n')
outFile.write('      seqtype = ' + args.seqtype + '\n')
outFile.write('    CodonFreq = ' + args.codonfreq + '\n')
outFile.write('\n')
outFile.write('        ndata = ' + args.ndata + '\n')
outFile.write('        clock = ' + args.clock + '\n')
outFile.write('       aaDist = ' + args.aaDist + '\n')
outFile.write('\n')
outFile.write('\n')
outFile.write('\n')
outFile.write('        model = ' + args.model + '\n')
outFile.write('\n')
outFile.write('\n')
outFile.write('\n')
outFile.write('\n')
outFile.write('\n')
outFile.write('\n')
outFile.write('      NSsites = ' + args.NSsites + '\n')
outFile.write('\n')
outFile.write('\n')
outFile.write('\n')
outFile.write('\n')
outFile.write('        icode = ' + args.icode + '\n')
outFile.write('        Mgene = ' + args.Mgene + '\n')
outFile.write('\n')
outFile.write('\n')
outFile.write('\n')
outFile.write('    fix_kappa = ' + args.fix_kappa + '\n')
outFile.write('        kappa = ' + args.kappa + '\n')
outFile.write('    fix_omega = ' + args.fix_omega + '\n')
outFile.write('        omega = ' + args.omega + '\n')
outFile.write('\n')
outFile.write('    fix_alpha = ' + args.fix_alpha + '\n')
outFile.write('        alpha = ' + args.alpha + '\n')
outFile.write('       Malpha = ' + args.Malpha + '\n')
outFile.write('        ncatG = ' + args.ncatG + '\n')
outFile.write('\n')
outFile.write('        getSE = ' + args.getSE + '\n')
outFile.write(' RateAncestor = ' + args.rateAncestor + '\n')
outFile.write('\n')
outFile.write('   Small_Diff = ' + args.Small_Diff + '\n')
outFile.write('    cleandata = ' + args.cleandata + '\n')
outFile.write('\n')
outFile.write('       method = ' + args.method + '\n')
outFile.close





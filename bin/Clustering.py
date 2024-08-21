#! /usr/bin/env python3
## Author: Sebastian Bechara
## Last Modified: 07.2024

import sys, os, argparse, shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from umap import UMAP

def argParse():
    '''
    Generates Parser using argparser
    Returns the arguments specified by command-line call and the parser itself
    '''
    des='\
This script will cluster the IRSs once per indivitual IES and once per KD experiment\n\
using different dimensionality reduction techniques (PCA, tSNE, UMAP).\n\
BE AWARE that dimensionality reduction with tSNE and UMAP is stochastic to a certain\n\
extent, which means that your results may differ slightly when repeated!!\n\
Per default it will only run a PCA on the KD level in 2 dimensions.\n\n\
Usage:\n\
      IRSclust.py [OPTIONS] -i <IRS_file> '
    parser = argparse.ArgumentParser(description=des,
                                    allow_abbrev=False,
                                    add_help=False,
                                    usage=argparse.SUPPRESS,
                                    formatter_class=argparse.RawTextHelpFormatter)
    req = parser.add_argument_group('Required')
    req.add_argument('-i','--IRS',
                     action='store',
                     metavar='',
                     type=str,
                     help='Specify the TSV-file containing the IRSs')

    optInOut = parser.add_argument_group('Optional - I/O options')
    optInOut.add_argument('-o','--outPath',
                          action='store',
                          metavar='',
                          type=str,
                          help='\
Specify the path used to save the generated plots. Will\n\
also be used for additional files. If it does not exist\n\
it will be created. Default is current directory')
    optInOut.add_argument('-k', '--KDs',
                          nargs='+',
                          metavar='',
                          help='\
Specify which KD (& ctrl) experiments to include in the\n\
analysis. If multiple, provide them seperated by whitespace.\n\
The spelling needs to be exactly as in the IRS file header.\n\
All are used per default')
    optInOut.add_argument('-f','--format',
                          action='store',
                          metavar='',
                          type=str,
                          help='\
Specify the format to save pictures in.\nAvailable formats: SVG [default], PNG, JPG')
    optInOut.add_argument('-D','--plot3D',
                          action='store_true',
                          help='Set if you wish a 3 dimensional plot. Default: False')
    optInOut.add_argument('-a', '--anim',
                          action='store_true',
                          help='\
Set this flag to generate an animated 3D scatterplot.\n\
It will rotate once around each axis. Overwrites and/or\n\
sets the following parameters: plot3D & dimensions. Default: False')

    optAdd = parser.add_argument_group('Optional - Additional tasks')
    optAdd.add_argument('-m','--method',
                        nargs='+',
                        metavar='',
                        help='\
Specify which method to Run. If multiple, provide them seperated\n\
by whitespace. Available methods: PCA, tSNE, UMAP [default]')
    optAdd.add_argument('-c','--clusterIES',
                        action='store_true',
                        help='\
Set this flag if you want to cluster the individual IESs in\n\
addition to the KD-experiments')
    optAdd.add_argument('-s', '--outScree',
                        action='store_true',
                        #metavar='',
                        #type=str,
                        help='Wether to create the scree plot in PCA. Default: False')
    optAdd.add_argument('-O','--outTSV',
                        action='store_true',
                        # metavar='',
                        # type=str,
                        help='Wether to write method coordinates to TSV. Default: False')
    optAdd.add_argument('-l','--loadings',
                        action='store_true',
                        help='Wether to add loadings of the Eigenvectors to the PCA biplot. Default: False')

    optParam = parser.add_argument_group('Optional - dimensionality reduction tuning')
    optParam.add_argument('-d','--dimensions',
                        #   nargs='+',
                          action='store',
                          metavar='',
                          type=int,
                          help='\
Specify an INT for the dimensions to calculate and to represent.\n\
For PCA this will set the number of the calculated principle component;\n\
only impacts the scree plot and if you want 3D biplots.\n\
For tSNE and UMAP this will be the dimensionality in which to represent\n\
your data; here only 2 or 3 makes sense. That means if >2 you will not get\n\
a 2D scatterplot and vice versa. (default: 2).')
    optParam.add_argument('-r','--randomState',
                          action='store',
                          metavar='',
                          type=int,
                          help='\
Specify an INT for reproducability. As dimensionality reduction\n\
has a stochastic element there is some random variance to be expected.\n\
When specifying a specific seed all subsequent iterations of the same\n\
data will produce the same results. Default: None')
    # optParam.add_argument('-M', '--metric',
    #                       action='store',
    #                       metavar='',
    #                       type=str,
    #                       help='Specify the metric to calculate distance in tSNE and UMAP. Available: euclidean, mahalanobis, chebyshev, ')
    optParam.add_argument('-p', '--perplexity',
                          action='store',
                          metavar='',
                          type=int,
                          help='\
Specify the perplexity used for tSNE. This is related to the nearest\n\
neighbours used in the algorithm.\nConcider a value between 5 and 50. Default is 30')
    optParam.add_argument('-n', '--neighbours',
                          action='store',
                          metavar='',
                          type=int,
                          help='\
Specify the number of neighbours to be used for UMAP.\n\
This is related to the balance between local and global structure\n\
detection (small - local; big - global).\nIs an INT between 2 and 100. Default: 15')
    optParam.add_argument('-M', '--min_dist',
                          action='store',
                          metavar='',
                          type=float,
                          help='\
Specify the minimum distance to space datapoints in the low\n\
dimensional space (UMAP only). This is related to the "packedness"\n\
of the resulting reduction (small - packed; big - spaced).\n\
Is an FLOAT between 0.001 and 0.99. Default: 0.1')
    

    opt = parser.add_argument_group('MISC')
    opt.add_argument('-h','--help',
                     action='help',
                     help='Print this message and exit\n\
For more info see documentation')
    opt.add_argument('-v', '--verbose',
                     action='store_true',
                     help='Set if you wish a verbose terminal output')

    return parser

def checkArgs(args, parser):
    ARGS={'IRS':None,
          'outPath':'.',
          'KDs':None,
          'F':'svg',
          'plot3D':False,
          'anim':False,
          'method':['PCA'],
          'cluster':False,
          'scree':False,
          'outTSV':False,
          'loadings':False,
          'dimensions':2,
          'seed':None,
          'perplexity':30,
          'neighbours':15,
          'min_dist':0.1,
          'V':False
          }
    if not args.IRS:
        print('ERROR!! No input file detected!')
        parser.print_help(sys.stderr)
        sys.exit()
    elif args.IRS.split('.')[-1].upper() != 'TSV':
        print('ERROR!! Wrong input format!')
        parser.print_help(sys.stderr)
        sys.exit()
    else: ARGS['IRS'] = args.IRS
    if args.outPath:
        if not os.path.isdir(args.outPath):
            os.mkdir(args.outPath)
            ARGS['outPath'] = args.outPath
        else: ARGS['outPath'] = args.outPath.rstrip('/') 
    if args.KDs:
        ARGS['KDs'] = args.KDs
    if args.format:
        if args.format.upper() not in ['SVG','PNG','JPG']:
            print('ERROR!! Format not supported!')
            parser.print_help(sys.stderr)
            sys.exit()
        else: ARGS['F'] = args.format.lower()
    if args.dimensions:
        if args.dimensions > 3:
            print('\tWARNING! If more than 3 dimensions are picked, no visual representation of your\ndimensionality reduction will be generated rather a TSV file continaing the values.')
            ARGS['dimensions'] = args.dimensions
            ARGS['outTSV'] = True
        else: ARGS['dimensions'] = args.dimensions
    if args.plot3D:
        ARGS['plot3D'] = True
        if not args.dimensions: 
            ARGS['dimensions'] = 3
        elif args.dimensions < 3: 
            ARGS['dimensions'] = 3
    if args.anim:
        ARGS['anim'] = True
        if not args.plot3D:
            ARGS['plot3D'] = True
        if args.dimensions:
            if not args.dimensions:
                ARGS['dimensions'] = 3
            elif args.dimensions < 3: 
                ARGS['dimensions'] = 3
        else: ARGS['dimensions'] = 3        
    if args.method:
        for meth in args.method:
            bad=[]
            if meth.lower() not in ['pca','tsne','umap']:
                bad.append(meth)
            if bad:
                print('ERROR!! Specified method not supported!')
                parser.print_help(sys.stderr)
                sys.exit()
        else: 
            ARGS['method']=args.method
    if args.clusterIES:
        ARGS['cluster']=True
    if args.outScree:
        ARGS['scree']=True
    if args.outTSV:
        ARGS['outTSV']=True
    if args.loadings:
        ARGS['loadings']=True

    if args.randomState:
        ARGS['seed']=args.randomState
    if args.perplexity:
        if not 5 <= args.perplexity <= 50:
            print('\tWARNING! Perplexity is outside this range {5, 50}.')
            ARGS['dimensions'] = args.dimensions
        else: ARGS['perplexity']=args.perplexity
    if args.neighbours:
        if not 2 <= args.neighbours <= 100:
            print('\tWARNING! Neighbours is outside this range {2, 100}.')
            ARGS['neighbours'] = args.neighbours
        else: ARGS['neighbours']=args.neighbours
    if args.min_dist:
        if not 2 <= args.min_dist <= 100:
            print('\tWARNING! Minimum distance for low dimensional embedding is outside this range {0.001, 0.99}.')
            ARGS['min_dist'] = args.min_dist
        else: ARGS['min_dist']=args.min_dist
    if args.verbose:
        ARGS['V'] = True
    return ARGS

def prepDF(IRSfile, sep='\t', KDs=None):
    '''
    Prepares the IRSs from the specified file for further processing.
    IES-ID is set as index and ID and Length are removed from the data.
    Returns two dataframes, one with the originial orientation and a 
    transposed one.
    '''
    df = pd.read_table(IRSfile, sep=sep, header = 0)
    df.index = df['ID']
    df = df.drop(['ID','Length'], axis = 1)
    if KDs:
        try:
            df = df[KDs]
        except:
            print('ERROR!! Specified KDs not in list! Maybe check the spelling?')
            parser.print_help(sys.stderr)
            sys.exit()
    df = df.fillna(0)
    featNames = df.columns
    df_T = df.transpose()
    return df, df_T, featNames

def scaleData(df, df_T):
    '''
    Scales the data using StrandardScaler() from sklearn. This essentially 
    Z-transforms the data, meaning each datapoint is substraced with the mean 
    and divided by the variance of the whole distribution.
    Returns two numpy arrays with the Z-transformed data (OG and transposed).  
    '''
    scaler = StandardScaler()
    X = scaler.fit_transform(df)
    scaler_T = StandardScaler()
    X_T = scaler_T.fit_transform(df_T)
    return X, X_T

def makeScatter(data, expVar=None, feature_names=None, 
    iesClust=False, loadings=False, method=None, out='.', F='svg', V=False):
    '''
    Function to generate the scatter plots for the dimensionality reductions.
    '''
    plt.figure(figsize=(12,9))
    
    if method=='PCA':
        method = 'PC'
    if iesClust:
        if V:
            print('\tCreating scatterplot for IESs')
        sns.scatterplot(data=data, x=f'{method} 1', y=f'{method} 2', color='darkgrey', alpha=0.05)
    else:
        if V:
            print('\tCreating scatterplot for KDs')
        sns.scatterplot(data=data, x=f'{method} 1', y=f'{method} 2', hue='labels')
    if loadings and method == 'PC':
        if iesClust:
            if V:
                print('\tAdding loadings')
            for ind, varName in enumerate(feature_names):
                plt.scatter( loadings[0][ind], loadings[1][ind])
                plt.arrow(0,0,
                            loadings[0][ind], loadings[1][ind],
                            color='black', width=0.0001, alpha=0.5)
                plt.text(loadings[0][ind], loadings[1][ind],
                        varName, ha='center', va='center')
        else:
            print('\tWARNING! Loadings will only be usefull for the individual IES clustering.\n\tSkipping . . .')

    elif loadings and method != 'PC':
        print('\tWARNING! Loadings were set but not a PCA.\n\tSkipping . . .')
    
    if method.startswith('PC'):
        plt.xlabel(f'PC1 (explained variance: {round(expVar[0]*100,2)} %)', fontsize=14)
        plt.ylabel(f'PC2 (explained variance: {round(expVar[1]*100,2)} %)', fontsize=14)
        # plt.title('PCA 1 & 2 Biplot', fontsize=16)
        if iesClust:
            plt.savefig(f'{out}/PCA_biPlot_IESs.{F}', bbox_inches='tight')
        else:
            plt.legend(ncols=7,bbox_to_anchor=(1.02,-0.06))
            plt.savefig(f'{out}/PCA_biPlot_KDs.{F}', bbox_inches='tight')
    else:
        plt.xlabel(f'{method} 1', fontsize=14)
        plt.ylabel(f'{method} 2', fontsize=14)
        # plt.title(f'{metho} 1 & 2 Plot', fontsize=16)
        if iesClust:
            plt.savefig(f'{out}/{method}_Plot_IESs.{F}', bbox_inches='tight')
        else:
            plt.legend(ncols=7,bbox_to_anchor=(1.02,-0.06))
            plt.savefig(f'{out}/{method}_Plot_KDs.{F}', bbox_inches='tight')

def scatter3D(X, iesClust=False, method=None, anim=False, outPath='.', F='svg', V=False): # X, 
    '''
    Create 3D scatterplot. Also makes it animated if wished
    '''
    fig = plt.figure(figsize=(9,9))
    ## for compatibility
    # ax = Axes3D(fig)
    # fig.add_axes(ax)
    ax = fig.add_subplot(111, projection='3d')

    if method == 'PCA':
        method = 'PC'
    
    if iesClust:
        if V:
            print('\tCreating 3D scatter plot for IESs')
        scat = ax.scatter(X[f'{method} 1'], X[f'{method} 2'], X[f'{method} 3'], s=0.5, linewidth=0)
    else:
        if V:
            print('\tCreating 3D scatter plot for KDs')
        scat = ax.scatter(X[f'{method} 1'], X[f'{method} 2'], X[f'{method} 3'], s=2)#, linewidth=0)

    ax.set_xlabel(f'{method} 1', fontsize=16)
    ax.set_ylabel(f'{method} 2', fontsize=16)
    ax.set_zlabel(f'{method} 3', fontsize=16)

    if anim:
        def update(frame):
            n_frames = 360 
            segment = n_frames // 4

            if frame < segment:
                # rotate z-axis
                angle = (frame  / segment) * 360
                ax.view_init(elev = 30, azim = angle)
            elif frame < 2 * segment:
                # rotate x-axis
                angle = ((frame-segment) / segment) * 360
                ax.view_init(elev = 30 + angle, azim = 360)
            elif frame < 3 * segment:
                # rotate y-axis
                angle = ((frame - 2 * segment ) / segment) * 360
                ax.view_init(elev = 390, azim = 360 + angle)
            else:
                # all rotations
                angle = ((frame - 3 * segment) / segment ) * 360
                ax.view_init(elev=390 + angle, azim = 720 + angle)
            return scat,

        anim = FuncAnimation(fig, update, frames=360, interval=100) 
        if iesClust:
            if V:
                print('\tCreating animated 3D scatterplot for IESs')
            anim.save(f'{outPath}/{["PCA" if method=="PC" else method][0]}_3D_plot_IES.gif', writer='pillow')
        else:
            if V:
                print('\tCreating animated 3D scatterplot for KDs')
            anim.save(f'{outPath}/{["PCA" if method=="PC" else method][0]}_3D_plot_KD.gif', writer='pillow')
    if iesClust:
        fig.savefig(f'{outPath}/{["PCA" if method=="PC" else method][0]}_3D_plot_IES.{F}')
    else:
        fig.savefig(f'{outPath}/{["PCA" if method=="PC" else method][0]}_3D_plot_KD.{F}')

def makeHex(data, method=None, out='.', F='svg'):
    '''
    Function to generate the hexbin plots for betrer visualisation.
    '''
    plt.figure(figsize=(12,9))
    plt.hexbin(data[f'{method} 1'], data[f'{method} 2'], gridsize=100, cmap='binary')
    cb = plt.colorbar(label='Counts')
    plt.title(f'{method} visualisation', fontsize=16)
    plt.xlabel(f'{method} 1',fontsize=14)
    plt.ylabel(f'{method} 2',fontsize=14)
    plt.savefig(f'{out}/{method}_hexbin.{F}', bbox_inches='tight')

def makeScree(expVar, outPath='.', iesClust=False, F='svg'):
    '''
    Creates Scree plots for PCA
    '''
    plt.figure(figsize=(9,6))
    plt.bar(range(1,len(expVar)+1), expVar)
    # plt.title(f'Scree plot')
    plt.xlabel('Principle Component')
    plt.ylabel('Explained variance')
    if iesClust:
        plt.savefig(f'{outPath}/PCA_Scree_IESlevel.{F}', bbox_inches='tight')
    else:
        plt.savefig(f'{outPath}/PCA_Scree.{F}', bbox_inches='tight')

def runPCA(X, X_T, components=None, outPath='.', format='svg', scree=False, 
           iesClust=False, loadings=False, feature_names=None, outTSV=False,
           plot3D=False, anim=False, V=False):
    '''
    Calculates PCA and generates biplots, screeplots and 3D plots
    '''
    if V:
        print('\tCalculating principle components on KD level')
    if components:
        pca_T=PCA(n_components=components)
    else:
        pca_T = PCA()
    pcaRes_T = pca_T.fit_transform(X_T)
    pcaPlot_T = pd.DataFrame(pcaRes_T)
    pcaPlot_T.columns = [f'PC {i}' for i in range(1,components+1)]#len(pcaPlot_T))]
    pcaPlot_T['labels'] = feature_names
    
    expVar_T = pca_T.explained_variance_ratio_
    if scree:
        if V:
            print('\tCreating Scree plot for KDs')
        makeScree(expVar_T, outPath, False, format)

    makeScatter(pcaPlot_T, expVar_T, feature_names, 
    False, loadings,'PCA', outPath, format, V)
    
    if plot3D:
        if components > 2:
            scatter3D(pcaPlot_T, False, 'PCA', anim, outPath, format, V)
        else:
            print('\tWARNING! 3D flag set, but less than 3 components specified.\n\tSkipping . . .')

    if outTSV:
        if V:
            print('\tWriting PCA results for KDs to file')
        pcaPlot_T.to_csv(f'{outPath}/PCA_{components}dim_vals_KDs.tsv', sep='\t', header=True, index=False)

    if iesClust:
        if V:
            print('\tCalculating principle components on IES level')
        if components:
            pca = PCA(n_components=components)
        else:
            pca = PCA()
        pcaRes = pca.fit_transform(X)
        pcaPlot = pd.DataFrame(pcaRes)
        pcaPlot.columns = [f'PC {i}' for i in range(1,components+1)]#len(pcaPlot))]
        
        expVar = pca.explained_variance_ratio_
        if scree:
            if V:
                print('\tCreating Scree plot for IESs')            
            makeScree(expVar, outPath, iesClust, format)
        
        makeScatter(pcaPlot, expVar, feature_names, 
        iesClust, loadings, 'PCA', outPath, format, V)
        
        if plot3D:
            if components > 2:
                scatter3D(pcaPlot, iesClust, 'PCA', anim, outPath, format, V)
            else:
                print('WARNING! 3D flag set, but less than 3 components specified.\nSkipping . . .')
        if outTSV:
            if V:
                print('\tWriting PCA results for IESs to file')
            pcaPlot.to_csv(f'{outPath}/PCA_{components}dim_vals_IESs.tsv', sep='\t', header=True, index=False)



def runTSNE(X, X_T, components=2, perplexity=30, iesClust=False, seed=None, outPath='.',
            F='svg', plot3D=False, anim=False, outTSV=False, feature_names=None, V=False):
    '''
    Calculates t-distributed neighbor embeding and creates plots
    '''
    if V:
        print(f'\tReducing dimensions to {components} components using tSNE on KD level')
    tsne_T = TSNE(n_components=components, 
                perplexity=perplexity, 
                random_state=seed)
    tsneRes_T = tsne_T.fit_transform(X_T)
    tsneDF_T = pd.DataFrame(tsneRes_T, columns=[f'tSNE {i+1}' for i in range(components)])
    tsneDF_T['labels'] = feature_names
    if components == 2: 
        makeScatter(tsneDF_T, None, None, False, False, 'tSNE', outPath, F, V)

    if plot3D:
        if components > 2:
            scatter3D(tsneDF_T, False, 'tSNE', anim, outPath, F, V)
        else:
            print('WARNING! 3D flag set, but less than 3 components specified.\nSkipping . . .')

    if outTSV:
        if V:
            print('\tWriting tSNE results for KDs to file')
        tsneDF_T.to_csv(f'{outPath}/tSNE_{components}dim_vals_KDs.tsv', sep='\t', header=True, index=False)
    
    if iesClust:
        if V:
            print(f'\tReducing dimensions to {components} components using tSNE on IES level')
        tsne = TSNE(n_components=components, 
                    perplexity=perplexity,  
                    random_state=None)
        tsneRes = tsne.fit_transform(X)
        tsneDF = pd.DataFrame(tsneRes, columns=[f'tSNE {i+1}' for i in range(components)])
        if components == 2:
            makeScatter(tsneDF, None, None, iesClust, False, 'tSNE', outPath, F, V)
            if V:
                print('\tCreating Hexbin plot')
            makeHex(tsneDF, 'tSNE', outPath, F)

        if plot3D:
            if components > 2:
                scatter3D(tsneDF, iesClust, 'tSNE', anim, outPath, F, V)
            else:
                print('WARNING! 3D flag set, but less than 3 components specified.\nSkipping . . .')

        if outTSV:
            if V:
                print('\tWriting tSNE results for IESs to file')
            tsneDF.to_csv(f'{outPath}/tSNE_{components}dim_vals_IESs.tsv', sep='\t', header=True, index=False)

def runUMAP(X, X_T, components=2, neighbours=15, minDist=0.1,
            iesClust=False, seed=None, feature_names=None, plot3D=False, anim=False, 
            outPath='.', F='svg', outTSV=False, V=False):
    '''
    Calculates uniform manifold approximation and projection, and creates plots
    '''
    if V:
        print(f'\tReducing dimensions to {components} components using UMAP on KD level')
    umap_T = UMAP(n_components=components, random_state=seed, n_neighbors=neighbours, min_dist=minDist)
    umapRes_T = umap_T.fit_transform(X_T)
    umapDF_T = pd.DataFrame(umapRes_T, columns=[f'UMAP {i+1}' for i in range(components)])
    umapDF_T['labels'] = feature_names
    
    if components == 2:
        makeScatter(umapDF_T, None, None, False, False, 'UMAP', outPath, F,V)

    if plot3D:
        if components > 2:
            scatter3D(umapDF_T, False, 'UMAP', anim, outPath, F, V)
        else:
            print('WARNING! 3D flag set, but less than 3 components specified.\nSkipping . . .')

    if outTSV:
        if V:
            print('\tWriting UMAP results for KDs to file')
        umapDF_T.to_csv(f'{outPath}/UMAP_{components}dim_vals_KDs.tsv', sep='\t', header=True, index=False)
    
    if iesClust:
        if V:
            print(f'\tReducing dimensions to {components} components using UMAP on IES level')        
        umap = UMAP(n_components=components, random_state=seed, n_neighbors=neighbours, min_dist=minDist)
        umapRes = umap.fit_transform(X)
        umapDF = pd.DataFrame(umapRes, columns=[f'UMAP {i+1}' for i in range(components)])
    
        if components == 2:
            makeScatter(umapDF, None, None, iesClust, False, 'UMAP', outPath, F, V)
            if V:
                print('\tCreating Hexbin plot')
            makeHex(umapDF, 'UMAP', outPath, F)

        if plot3D:
            if components > 2:
                scatter3D(umapDF, iesClust, 'UMAP', anim, outPath, F, V)
            else:
                print('WARNING! 3D flag set, but less than 3 components specified.\nSkipping . . .')

        if outTSV:
            if V:
                print('\tWriting UMAP results for IESs to file')
            umapDF.to_csv(f'{outPath}/UMAP_{components}dim_vals_IESs.tsv', sep='\t', header=True, index=False)

def getTime():
    return str(datetime.now()).split('.')[0].replace('-','.').replace(' ','|')

def main(args=None):
    parser = argParse()
    num=1
    if args:
        num=2
    if not args:
        args = parser.parse_args()
    if len(sys.argv) == num:
        parser.print_help(sys.stderr)
        sys.exit()
    
    ARGS = checkArgs(args, parser)
    if ARGS['V']:
        print(f'\n\n[{getTime()}] - Starting clustering analysis')
        
        print(f'\
\tThe following Parameters were specified:\n\
\t\tIRS file: {ARGS["IRS"]}\n\
\t\tOutput directory: {ARGS["outPath"]}\n\
\t\tMethods used for clustering: {", ".join(ARGS["method"])}\n\
\t\tClustering on IES level: {ARGS["cluster"]}\n\
\t\tKDs to consider: {[", ".join(ARGS["KDs"]) if ARGS["KDs"] else "None - using all"][0]}\n\
\t\tPicture format: {ARGS["F"]}\n\
\t\t3D plots: {ARGS["plot3D"]}\n\
\t\tAnimated 3D plots: {ARGS["anim"]}\n\
\t\tOutput reduced values: {ARGS["outTSV"]}\n\
\t\tThe folowing parameters are directly impacting the algorithms:\n\
\t\t[All] Dimensions to reduce to: {ARGS["dimensions"]}\n\
\t\t[PCA] Generate Scree plot: {ARGS["scree"]}\n\
\t\t[PCA] Add loadings: {ARGS["loadings"]}\n\
\t\t[tSNE & UMAP] RNG seed: {ARGS["seed"]}\n\
\t\t[tSNE] Perplexity: {ARGS["perplexity"]}\n\
\t\t[UMAP] Neighbours: {ARGS["neighbours"]}\n\
\t\t[UMAP] Minimum distance: {ARGS["min_dist"]}')        
            
    else:
        print(f'\n\n[{getTime()}] - Parsing and checking arguments')
        # ARGS = checkArgs(args, parser)
    print(f'[{getTime()}] - Setting up data')
    if ARGS['V']:
        print('\tReading IRS-file')
    df, df_T, featNames = prepDF(ARGS['IRS'], sep='\t', KDs=ARGS['KDs'])
    if ARGS['V']:
        print('\tScaling Data')
    X, X_T = scaleData(df, df_T)

    for method in ARGS['method']:
        if method.upper() == 'PCA':
            print(f'[{getTime()}] - Running PCA')
            runPCA(X, X_T, components=ARGS['dimensions'], outPath=ARGS['outPath'], 
            format=ARGS['F'], scree=ARGS['scree'], iesClust=ARGS['cluster'], 
            loadings=ARGS['loadings'], feature_names=featNames,
            outTSV=ARGS['outTSV'], plot3D=ARGS['plot3D'], anim=ARGS['anim'], V=ARGS['V'])

        if method.upper() == 'TSNE':
            print(f'[{getTime()}] - Running tSNE')
            runTSNE(X, X_T, components=ARGS['dimensions'], perplexity=ARGS['perplexity'],
            iesClust=ARGS['cluster'], seed=ARGS['seed'], outPath=ARGS['outPath'],
            F=ARGS['F'], plot3D=ARGS['plot3D'], anim=ARGS['anim'], outTSV=ARGS['outTSV'],
            feature_names=featNames, V=ARGS['V'])

        if method.upper() == 'UMAP':
            print(f'[{getTime()}] - Running UMAP')
            runUMAP(X, X_T, components=ARGS['dimensions'], neighbours=ARGS['neighbours'],
            minDist=ARGS['min_dist'], iesClust=ARGS['cluster'], seed=ARGS['seed'], 
            feature_names=featNames, plot3D=ARGS['plot3D'], anim=ARGS['anim'], 
            outPath=ARGS['outPath'], F=ARGS['F'], outTSV=ARGS['outTSV'], V=ARGS['V'])
    if ARGS['V']:
        if ARGS['outPath'] == '.':
            print(f'[{getTime()}] - Done!!\n\nYou can find your results in this directory:\n\t{ARGS["outPath"]} - meaning the current working directory')
            print('\tGenerated Files:')
            for file in os.listdir(ARGS['outPath']):
                print(f'\t\t - {file}')
        else:
            print(f'[{getTime()}] - Done!!\n\nYou can find your results in this directory:\n\t{ARGS["outPath"]}')
            print('\tGenerated Files:')
            for file in os.listdir(ARGS['outPath']):
                print(f'\t\t - {file}')
    else:
        print(f'[{getTime()}] - Done!!')

if __name__ == '__main__':
    main()

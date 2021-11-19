import os
import sys
from matplotlib import font_manager
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
import mpl_scatter_density as mplSD
# colMap = cm.get_cmap('Set1', 5)
# print(colMap.colors)
# for color in colMap.colors:
#     print(mpl.colors.to_hex(color))

#colorblind palette
dColors = {
    'MCF102A': '#029E73',
    'SKBR3': '#D55E00',
    'MCF7': '#CC78BC',
    'MDAMB231': '#BAB029',
    'MDAMB361': '#0173B2'
}

# plt.style.use('engro2')


def bareStyle(figure=None, axes=None):
    # axes.grid(False)
    axes.grid(True)
    # axes.set_frame_on(False)
    axes.set_frame_on(True)
    axes.set_aspect('auto')
    axes.tick_params(axis='both',
                     which='both',
                     bottom=False,
                     top=False,
                     left=False,
                     right=False,
                     labelbottom=False,
                     labeltop=False,
                     labelleft=False,
                     labelright=False)
                     # grid_color = 'white',
                     # grid_alpha = 1)


def tSneTagLabels(axes, fileName, dizTagColors,):
    dfTags = pd.read_csv(fileName, sep='\t')
    # print(dfTags)
    for tag in dfTags.itertuples(name='row'):
        axes.text(x=tag.xTag,
                  y=tag.yTag,
                  s=tag.Tag,
                  color=dizTagColors[tag.Tag],
                  transform=axes.transAxes,
                  fontsize=18)


def scatterOnTsne(df,
                  fluxName,
                  colorMap='viridis',
                  sTitle='',
                  tagFileName='',
                  dizTagColors=dColors,
                  file2save='sctrOnTsne.pdf'):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_axes([0.00, 0.05, 1.0, 0.91])
    # ax = fig.add_axes([0.0, 0.0, 0.9, .95])
    scPlt = ax.scatter(df['tsne-2d-one'],
                       df['tsne-2d-two'],
                       c=df[fluxName],
                       cmap=colorMap,
                       s=4,
                       alpha=0.3,
                       edgecolors='none')
    cbar = fig.colorbar(scPlt, drawedges=True, shrink=0.5)
    cbar.dividers.set_visible(False)
    cbar.outline.set_visible(False)
    # cbar.solids.set_edgecolor("face")
    # cbar = fig.colorbar(scPlt, shrink=0.8)
    bareStyle(axes=ax)
    tSneTagLabels(ax, tagFileName, dizTagColors)
    ax.set_title(sTitle)
    plt.savefig(file2save)
    plt.close()


def tagsOnTsne(df,
               tagColumn,
               dizTagColors=dColors,
               tagFileName='',
               sTitle='',
               file2save='tSNEtags.pdf'):
    fig = plt.figure(figsize=(10, 10))

    # ax = fig.add_axes([0.05, 0.05, 0.95, 0.95])
    ax = fig.add_axes([0.05, 0.05, 0.8, 0.8])

    ax.set_axisbelow(True)

    plt.grid(color='w', linestyle='solid')

    # hide axis spines
    for spine in ax.spines.values():
        spine.set_visible(False)

    for tagColor in dizTagColors:
        df2Plot = df[df[tagColumn] == tagColor]
        ax.scatter(df2Plot['tsne-2d-one'],
                   df2Plot['tsne-2d-two'],
                   label=tagColor,
                   c=dizTagColors[tagColor],
                   s=4,
                   alpha=0.3,
                   edgecolors='none')

    bareStyle(axes=ax)
    ax.set_facecolor('#eaeaf2')
    # ax.xaxis.set_visible(False)
    # ax.yaxis.set_visible(False)
    tSneTagLabels(ax, tagFileName, dizTagColors)

    ax.set_xlabel('tSNE_1', fontsize = 20, color = 'grey')
    ax.set_ylabel('tSNE_2', fontsize = 20, color = 'grey')

    plt.title(sTitle)
    # plt.legend()
    plt.savefig(file2save)
    plt.close()

def tagsOnTsne_plotMultiple(df, df2,
               tagColumn,
               dizTagColors=dColors,
               tagFileName='',
               sTitle='',
               file2save='tSNEtags.pdf'):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_axes([0.05, 0.05, 0.95, 0.95])
    for tagColor in dizTagColors:
        df2Plot = df[df[tagColumn] == tagColor]
        ax.scatter(df2Plot['tsne-2d-one'],
                   df2Plot['tsne-2d-two'],
                   label=tagColor,
                   c=dizTagColors[tagColor],
                   s=4,
                   alpha=0.3,
                   edgecolors='none')
        df2Plot2 = df2[df2[tagColumn] == tagColor]
        ax.scatter(df2Plot2['tsne-2d-one'],
                   df2Plot2['tsne-2d-two'],
                   label=tagColor,
                   c=dizTagColors[tagColor],
                   s=4,
                   alpha=0.3,
                   edgecolors='none')
        # df3Plot3 = df3[df3[tagColumn] == tagColor]
        # ax.scatter(df3Plot3['tsne-2d-one'],
        #            df3Plot3['tsne-2d-two'],
        #            label=tagColor,
        #            c=dizTagColors[tagColor],
        #            s=4,
        #            alpha=0.3,
        #            edgecolors='none')
        # df4Plot4 = df4[df4[tagColumn] == tagColor]
        # ax.scatter(df4Plot4['tsne-2d-one'],
        #            df4Plot4['tsne-2d-two'],
        #            label=tagColor,
        #            c=dizTagColors[tagColor],
        #            s=4,
        #            alpha=0.3,
        #            edgecolors='none')
    bareStyle(axes=ax)
    tSneTagLabels(ax, tagFileName, dizTagColors)
    plt.title(sTitle)
    # plt.legend()
    plt.savefig(file2save)
    plt.close()


def histoDistribs(df,
                  fluxName,
                  binDistrib,
                  tagColumn,
                  dizTagColors,
                  sTitle='',
                  file2save='histioDistribs.pdf'):
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111)
    for tag in dizTagColors:
        ax.hist(df[df[tagColumn] == tag][fluxName],
                bins=binDistrib,
                alpha=0.2,
                label=tag,
                density=False,
                color=dColors[tag])
    ax.set_xlabel('Flux value')
    ax.set_ylabel('Counts')
    plt.title(sTitle)
    plt.savefig(file2save)
    plt.close()


def scatAndHistOnTsne(dfwTsne, df,
                      fluxName,
                      colorMap='viridis',
                      sTitle='',
                      tagFileName='',
                      tagColumn='',
                      dizTagColors=dColors,
                      binDistrib='',
                      file2save='sctrOnTsne.pdf'):
    # fig = plt.figure(figsize=(10, 16))
    # gs = gridspec.GridSpec(nrows=16, ncols=1, figure=fig, wspace=0, hspace=0)
    # axScat = fig.add_subplot(gs[:11, 0])
    # axHist = fig.add_subplot(gs[11:14, 0], frameon=False)
    # axCBar = fig.add_subplot(gs[-1, 0], sharex=axHist, frameon=False)

    fig = plt.figure(figsize=(20, 10))
    gs = gridspec.GridSpec(nrows=10, ncols=20, figure=fig, wspace=0, hspace=0)
    # axScat = fig.add_subplot(gs[:, :10])
    axScat = fig.add_subplot(gs[:, :8])
    axHist = fig.add_subplot(gs[:8, 11:], frameon=False)
    axCBar = fig.add_subplot(gs[9, 11:], sharex=axHist, frameon=False)

    font_dirs = ['/home/mdifilippo/Roboto']
    font_files = font_manager.findSystemFonts(fontpaths=font_dirs)

    for font_file in font_files:
        font_manager.fontManager.addfont(font_file)

    mpl.rcParams['font.sans-serif'] = "Roboto"
    mpl.rcParams['font.family'] = "Roboto"

    # ax = fig.add_axes([0.0, 0.0, 0.9, .95])
    scPlt = axScat.scatter(dfwTsne['tsne-2d-one'],
                           dfwTsne['tsne-2d-two'],
                           c=dfwTsne[fluxName],
                           cmap=colorMap,
                           s=4,
                           alpha=0.7
                           # , edgecolors='none'
                           )
    # axHist = fig.add_subplot(3, 1, 2, frameon=False)
    for tag in dizTagColors:
        axHist.hist(df[df[tagColumn] == tag][fluxName],
                    bins=binDistrib,
                    alpha=0.2,
                    label=tag,
                    density=False,
                    color=dColors[tag])
        axHist.set_xlabel('Flux value', fontsize = 22)
        axHist.set_ylabel('Counts', fontsize = 22)

        axHist.tick_params(axis='both', labelsize=22)
        # axHist.set_yticklabels(fontsize = 22)
        # axHist.set_xticklabels(fontsize = 22)
    plt.xticks(fontsize= 22)
    plt.yticks(fontsize= 22)
    # axCBar = fig.add_subplot(3, 1, 3, sharex=axHist)
    cbar = fig.colorbar(scPlt,
                        cax=axCBar,
                        orientation='horizontal',
                        drawedges=True,
                        fraction=0.1)
    cbar.dividers.set_visible(False)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(labelsize=22)
    # cbar.solids.set_edgecolor("face")
    # cbar = fig.colorbar(scPlt, shrink=0.8)
    bareStyle(axes=axScat)
    tSneTagLabels(axScat, tagFileName, dizTagColors)
    # axScat.set_title(sTitle, loc='right', fontdict={'color': '#808080'}, fontsize = 25)
    axHist.set_title(sTitle, loc='left', fontdict={'color': '#808080'}, fontsize = 35)
    fig.subplots_adjust(hspace=0)
    plt.savefig(file2save)
    plt.close()


def restrict_q5_q95(df, fluxName):
    (q5, q95) = df[fluxName].quantile([0.05, 0.95])
    # print(q5, q95)
    inRange = (df[fluxName] >= q5) & (df[fluxName] <= q95)
    # print('n Points in range', inRange[inRange == True].count(),
    #       'n Points left Out', inRange[inRange == False].count())
    dfOut = df[inRange]
    nBins = np.ceil(2 * dfOut.shape[0]**(1 / 3))
    bins = np.linspace(q5, q95, num=int(nBins + 1))
    # print(nBins, bins)
    return (dfOut, nBins, bins)


def separationScoreDF(dfFluxes, fluxName, lcellLines, histoBins):
    dfScores = pd.DataFrame()
    dfHistos = pd.DataFrame()
    dCL2Subtract = {}
    for line in lcellLines:
        # print(dfFluxes['cellLine'] == line)
        sCLine = dfFluxes[dfFluxes['cellLine'] == line][fluxName]
        histoCellLine, binEdges = np.histogram(sCLine, histoBins)
        dfHistos[line] = histoCellLine
        lCtmp = lcellLines.copy()
        lCtmp.remove(line)
        dCL2Subtract[line] = lCtmp
        # print(histoCellLine)
    # print(dfHistos)
    # print(dCL2Subtract)
    for line in lcellLines:
        scoreColName = line + '_score'
        dfHistos[scoreColName] = dfHistos[line] - \
            dfHistos[dCL2Subtract[line]].sum(axis=1)
        dfHistos[scoreColName] /= dfHistos[line].sum()
    dfHistos.where(dfHistos > 0, other=0, inplace=True)
    dfScores[fluxName] = dfHistos.iloc[:, 5:].sum()
    totalScore = dfHistos.iloc[:, 5:].sum().sum()
    dfScores.loc['totalScore', fluxName] = totalScore
    return dfScores


def biomassVSexchange(dfFluxes,
                      bioMFName,
                      ExRxnName,
                      lcellLines,
                      dCellLinesColors=dColors,
                      sTitle='',
                      file2save=''):
    fig = plt.figure(figsize=(10, 10))
    # ax = fig.add_axes([0.00, 0.05, 1.0, 0.91])
    ax = fig.add_axes([0.17, 0.1, 0.8, .85])
    exMax = dfFluxes[ExRxnName].min()
    exMax *= -1.0
    ax.fill_between([0, exMax], [0, exMax], fc='#C0C0C0', alpha=0.2)
    for cLine in lcellLines:
        dfExs = -1.0 * dfFluxes[dfFluxes.cellLine == cLine][ExRxnName]
        ax.scatter(dfExs,
                   dfFluxes[dfFluxes.cellLine == cLine][bioMFName],
                   c=dCellLinesColors[cLine],
                   s=4,
                   alpha=0.5,
                   edgecolors='none')
    ax.set_title(sTitle)
    ax.set_xlabel('Exchange Rxn value')
    ax.set_ylabel('Biomass Flux value')
    plt.savefig(file2save)
    plt.close()

    # lMethods = ['sqrt', 'sturges', 'rice', 'scott', 'fd', 'doane']
    # for m in lMethods:
    #     fig = plt.figure(figsize=(13.3, 10))
    #     ax = fig.add_subplot(111)
    #     for cellLine in lcellLines:
    #         ax.hist(dfToPlot[dfToPlot.cellLine == cellLine]
    #                 [flux], bins=m, alpha=0.3, label=cellLine,
    #                 color=dvl.dColors[cellLine])
    #     ax.set_xlabel('Flux values')
    #     ax.set_ylabel('Frequency')
    #     plt.title('random Sampling - ' + m)
    #     plt.savefig(os.path.join(
    #         FIGUREDIR, basenameRS + '_hist_' + m + '_all.pdf'))
    #     plt.close()


def biomassVSexchangeNeg(dfFluxes,
                         bioMFName,
                         ExRxnName,
                         coeffExInBio,
                         lcellLines,
                         dCellLinesColors=dColors,
                         sTitle='',
                         file2save=''):
    fig = plt.figure(figsize=(10, 10))
    # ax = fig.add_axes([0.00, 0.05, 1.0, 0.91])
    ax = fig.add_axes([0.17, 0.1, 0.8, .85])
    exMaxVal = dfFluxes[ExRxnName].max()
    exMaxIdx = dfFluxes[ExRxnName].idxmax()
    exMinVal = dfFluxes[ExRxnName].min()
    exMinIdx = dfFluxes[ExRxnName].idxmin()
    #  biomMax = dfFluxes.loc[exMaxIdx, bioMFName]
    #  biomMin = dfFluxes.loc[exMinIdx, bioMFName]
    biomMax = coeffExInBio * -1.0 * exMaxVal
    biomMin = coeffExInBio * -1.0 * exMinVal
    print('EX_rxn', ExRxnName, 'V_ExInBiom', bioMFName, 'coeffExInBio',
          coeffExInBio)
    print('idx', 'V_ex', 'V_biom')
    print('Max', exMaxIdx, exMaxVal, biomMax)
    print('Min', exMinIdx, exMinVal, biomMin)
    if exMinVal >= 0:
        biomMax *= -1.0
        ax.fill_between([0, exMaxVal], [0, biomMax], fc='#C0C0C0', alpha=0.2)
    else:
        ax.fill_between([exMinVal, 0], [biomMin, 0], fc='#C0C0C0', alpha=0.2)
    for cLine in lcellLines:
        dfExs = dfFluxes[dfFluxes.cellLine == cLine][ExRxnName]
        ax.scatter(dfExs,
                   dfFluxes[dfFluxes.cellLine == cLine][bioMFName],
                   c=dCellLinesColors[cLine],
                   s=4,
                   alpha=0.5,
                   edgecolors='none',
                   label=cLine)
    ax.set_title(sTitle)
    ax.set_xlabel('Exchange Rxn value')
    ax.set_ylabel('Biomass Flux value')
    plt.legend()
    plt.savefig(file2save)
    plt.close()


def vExcngInBiomass(dfFluxes,
                    bioMFName,
                    ExRxnName,
                    coeffExInBio,
                    lcellLines,
                    dCellLinesColors=dColors,
                    sTitle='',
                    file2save=''):
    fig = plt.figure(figsize=(10, 10))

    font_dirs = ['/home/mdifilippo/Roboto']
    font_files = font_manager.findSystemFonts(fontpaths=font_dirs)

    for font_file in font_files:
        font_manager.fontManager.addfont(font_file)

    mpl.rcParams['font.sans-serif'] = "Roboto"
    mpl.rcParams['font.family'] = "Roboto"

    # ax = fig.add_axes([0.00, 0.05, 1.0, 0.91])
    # ax = fig.add_axes([0.17, 0.1, 0.8, .85])
    ax = fig.add_axes([0.1, 0.1, 0.8, .85])
    exMaxVal = dfFluxes[ExRxnName].max()
    exMaxIdx = dfFluxes[ExRxnName].idxmax()
    exMinVal = dfFluxes[ExRxnName].min()
    exMinIdx = dfFluxes[ExRxnName].idxmin()
    #  biomMax = dfFluxes.loc[exMaxIdx, bioMFName]
    #  biomMin = dfFluxes.loc[exMinIdx, bioMFName]i
    biomMax = dfFluxes[bioMFName].max()
    biomMin = dfFluxes[bioMFName].min()
    print('EX_rxn', ExRxnName, 'V_ExInBiom', bioMFName, 'coeffExInBio',
          coeffExInBio)
    print('idx', 'V_ex', 'V_biom')
    print('Max', exMaxIdx, exMaxVal, biomMax)
    print('Min', exMinIdx, exMinVal, biomMin)
    if exMinVal >= 0:
        ax.fill_betweenx([0, biomMax], [0, biomMax], fc='#C0C0C0', alpha=0.2)
        ax.plot([0, biomMax], [0, biomMax], lw=1.0, c='#C0C0C0', alpha=0.8)
    else:
        ax.fill_betweenx([biomMax, 0], [-biomMax, 0], fc='#C0C0C0', alpha=0.2)
        ax.plot([-biomMax, 0], [biomMax, 0], lw=1.0, c='#C0C0C0', alpha=0.8)
    for cLine in lcellLines:
        dfExs = dfFluxes[dfFluxes.cellLine == cLine][ExRxnName]
        ax.scatter(dfExs,
                   dfFluxes[dfFluxes.cellLine == cLine][bioMFName],
                   c=dCellLinesColors[cLine],
                   s=4,
                   alpha=0.5,
                   edgecolors='none',
                   label=cLine)
    ax.set_title(sTitle, fontsize = 30)
    ax.set_xlabel('Exchange Rxn value', fontsize = 22)
    ax.set_ylabel('Biomass Flux value', fontsize = 22)
    ax.tick_params(axis='both', labelsize=22)
    ax.legend(labelcolor=dColors.values(), scatterpoints=0, frameon=False, fontsize = 22)

    lxticks = plt.xticks()[0]
    lNewValues_x = [int(lxticks[0])]
    for r in lxticks[1:]:
        lNewValues_x.append(round(r,3))
    plt.xticks(lNewValues_x,lNewValues_x, fontsize= 22)

    lyticks = plt.yticks()[0]
    lNewValues_y = [int(lyticks[0])]
    for r in lyticks[1:]:
        lNewValues_y.append(round(r,4))
    plt.yticks(lNewValues_y,lNewValues_y, fontsize= 22)

    plt.savefig(file2save)
    plt.close()


def vExcngInBiomassScaterDensity(dfFluxes,
                    bioMFName,
                    ExRxnName,
                    coeffExInBio,
                    lcellLines,
                    dCellLinesColors=dColors,
                    sTitle='',
                    file2save=''):
    fig = plt.figure(figsize=(10, 10))
    # ax = fig.add_axes([0.00, 0.05, 1.0, 0.91])
    ax = fig.add_axes([0.17, 0.1, 0.8, .85], projection='scatter_density')
    exMaxVal = dfFluxes[ExRxnName].max()
    exMaxIdx = dfFluxes[ExRxnName].idxmax()
    exMinVal = dfFluxes[ExRxnName].min()
    exMinIdx = dfFluxes[ExRxnName].idxmin()
    #  biomMax = dfFluxes.loc[exMaxIdx, bioMFName]
    #  biomMin = dfFluxes.loc[exMinIdx, bioMFName]i
    biomMax = dfFluxes[bioMFName].max()
    biomMin = dfFluxes[bioMFName].min()
    print('EX_rxn', ExRxnName, 'V_ExInBiom', bioMFName, 'coeffExInBio',
          coeffExInBio)
    print('idx', 'V_ex', 'V_biom')
    print('Max', exMaxIdx, exMaxVal, biomMax)
    print('Min', exMinIdx, exMinVal, biomMin)
    if exMinVal >= 0:
        ax.fill_betweenx([0, biomMax], [0, biomMax], fc='#C0C0C0', alpha=0.2)
        ax.plot([0, biomMax], [0, biomMax], lw=1.0, c='#C0C0C0', alpha=0.8)
    else:
        ax.fill_betweenx([biomMax, 0], [-biomMax, 0], fc='#C0C0C0', alpha=0.2)
        ax.plot([-biomMax, 0], [biomMax, 0], lw=1.0, c='#C0C0C0', alpha=0.8)
    for cLine in lcellLines:
        dfExs = dfFluxes[dfFluxes.cellLine == cLine][ExRxnName]
        ax.scatter_density(dfExs,
                   dfFluxes[dfFluxes.cellLine == cLine][bioMFName],
                   color=dCellLinesColors[cLine],
                   s=4,
                   alpha=0.5,
                   edgecolors='none',
                   label=cLine)
    ax.set_title(sTitle)
    ax.set_xlabel('Exchange Rxn value')
    ax.set_ylabel('{} contribution to biomass Flux value'.format(sTitle))
    ax.legend(labelcolor=dColors.values(), scatterpoints=0, frameon=False)
    plt.savefig(file2save)
    plt.close()

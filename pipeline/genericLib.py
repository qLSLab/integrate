import os
import sys
import pandas as pd
import numpy as np
import re
import time


def setWorkingDirs(dataDir='rawData', outDir=None):
    """Set the working directories dataDir, outDir, mapsDir, reportDir.

    If arguments are omitted current working directory replace them.

    Keyword arguments:
     dataDir: The directory to look for the input files.
     outDir:  The directory to write the output files.

    Return:
    dataDir, outDir, reportDir, modelDir, logDir, figureDir, scriptsDir, mapDir
    """
    cwDir = os.getcwd()
    if dataDir is None:
        dataDir = cwDir + os.sep
    else:
        if not os.path.exists(dataDir):
            # print'the input directory', dataDir, 'does not exists. Fix it.'
            sys.exit()
        else:
            dataDir = os.path.abspath(dataDir)
            if dataDir[-1] != os.sep:
                dataDir += os.sep
    if outDir is None:
        outDir = cwDir
    else:
        outDir = os.path.abspath(outDir) + os.sep
    reportDir = pathJoinCheck('outputs', outDir)
    modelDir = pathJoinCheck('models', outDir)
    figureDir = pathJoinCheck('figures', outDir)
    pipelineDir = pathJoinCheck('pipeline', outDir)
    return dataDir, outDir, reportDir, modelDir, figureDir, pipelineDir

def pathJoinCheck(dir2add, rootPath='.'):
    """Check the existence of a path and build it if necessary."""
    path = os.path.join(rootPath, dir2add)
    if not os.path.exists(path):
        os.makedirs(path)
    return path

def openFile(name, estensione, dir=None):
    cwDir = os.getcwd()
    if dir is None:
        dir = cwDir + os.sep
    fileName = os.path.join(dir, name + '.' + estensione)
    # printfileName
    stream = open(fileName, mode='w')
    return stream

def writeLineByLineToFile(stream, dataToWrite, spaziatore):
    stream.write(spaziatore.join(str(s) for s in dataToWrite) + '\n')

def extractRegexFromItem(item, regex):
    sItem = pd.Series(data=item)
    if sItem.empty == False:
        dfItem = sItem.str.extractall(regex)
    else:
        dfItem = []
    return dfItem


def intersect(a, b):
	""" return the intersection of two lists """
	return list(set(a) & set(b))


##########################
def getEuclideanDistance(v1, v2):
    return np.round(np.sqrt(np.sum((v1 - v2) ** 2)), decimals = 5)


def reduceValueByGivenPercentage(value, percentage):
    variation = (value * percentage) / 100
    return value - variation

def setPath(path):
    """Check the existence of a path and build it if necessary."""
    if not os.path.exists(path):
        os.makedirs(path)


def createEmptyDictWithKeyTarget(keysTarget):
    dEmpty = {}
    for k in keysTarget:
        dEmpty[k] = None
    return dEmpty

def unique(a):
	""" return the list with duplicate elements removed """
	return list(set(a))

def union(a, b):
	""" return the union of two lists """
	return list(set(a) | set(b))

def difference(a, b):
    """ return the difference of two lists """
    return list(set(a) - set(b))


def splitStringEveryNChars(string, n):
    return [string[i:i+n] for i in range(0, len(string))]


def getTimeStamp():
    return time.strftime('%Y%m%d%H%M%S', time.localtime())


def extractPercentage(element):
    el_str = element.split(' ')[0].strip()
    dfElement = extractRegexFromItem(el_str, r"([0-9]+)")
    ldfElement = list(dfElement[0])
    return int(ldfElement[0])

def logFileOpen(logDIR=None, timeStamp=None, aim = None):
    cwDir = os.getcwd()
    if logDIR is None:
        logDIR = cwDir + os.sep
    if timeStamp is None:
        timeStamp = getTimeStamp()
    logFileName = os.path.join(logDIR, timeStamp + '_' + aim + '.log')
    logStream = open(logFileName, mode='w')
    return logStream

def toLog(logStream, string):
    logStream.write(string + '\n')


# 1 a + 2 b -> 3-c
def dizReaProd(s):
    termini = s.strip().split('+')
    # print('termini', termini
    diz = {}
    for termine in termini:
        coeffMetabo = termine.split()
        # print('coeffMetabo', coeffMetabo
        coeff = coeffMetabo[0]
        if isCoeff(coeff) is True:
            coeff = float(coeffMetabo[0])
            metabolita = ' '.join(coeffMetabo[1:])
        else:
            metabolita = ' '.join(coeffMetabo)
            coeff = 1.0
        # print('coeff', coeff
        # print('metabolita', metabolita
        diz[metabolita] = coeff
    return diz

def isCoeff(s):
    """Determine if a string splitted on the spaces the first element is the
    stoichiometric coefficient or not.
    Example: if string is "2 A" return True; if string is "A" it returns False"""
    answer = re.match('((\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)$', s)
    if answer is not None:
        # print('\t match ', answer.group(0)
        return True
    else:
        return False

# def findAndReplaceStr(stringa, substringa, newStringa):
#     pos = stringa.find(substringa)
#     print('pos\t', pos)
#     badPositions = []
#     while pos != -1 and pos not in badPositions:
#         if pos > 0 and (stringa[pos - 1].isdigit() == False and stringa[pos - 1].isalpha() == False and stringa[pos - 1] != '-') \
#                 and pos + len(substringa) < len(stringa) and (stringa[pos + len(substringa)].isdigit() == False and stringa[pos + len(substringa)].isalpha() == False and stringa[pos + len(substringa)] != '-'):
#             stringa = stringa[:pos] + str(newStringa) + stringa[pos + len(substringa):]
#         elif pos == 0 and pos + len(substringa) < len(stringa) and (stringa[pos + len(substringa)].isdigit() == False and \
#             stringa[pos + len(substringa)].isalpha() == False and stringa[pos + len(substringa)] != '-'):
#             stringa = stringa[:pos] + str(newStringa) + stringa[pos + len(substringa):]
#         elif pos > 0 and pos + len(substringa) == len(stringa) and (stringa[pos - 1].isdigit() == False and \
#             stringa[pos - 1].isalpha() == False and stringa[pos - 1] != '-'):
#             stringa = stringa[:pos] + str(newStringa) + stringa[pos + len(substringa):]
#         elif pos == 0 and pos + len(substringa) == len(stringa):
#             stringa = str(newStringa)
#         else:
#             badPositions.append(pos)
#         pos = stringa.find(substringa)
#         print('new pos\t', pos)
#         print('badPositions\t', badPositions,'\n', stringa, '\n')
#     return stringa

def findAndReplaceStr(stringa, substringa, newStringa):
    # pos = stringa.find(substringa)
    lpos = [i.start() for i in re.finditer(substringa, stringa)]
    # print('iniziale Stringa\t', stringa)
    # print('lpos\t', lpos)
    badPositions = 0
    lenNewStringa = len(newStringa)
    # print('lenNewStringa\t', lenNewStringa)
    k = 0
    # print('k iniziale\t', k)
    while len(lpos) != badPositions and k < len(lpos):
        # print('current element\t', lpos[k])
        if lpos[k] > 0 and (stringa[lpos[k] - 1].isdigit() == False and stringa[lpos[k] - 1].isalpha() == False and stringa[lpos[k] - 1] != '-') \
                and lpos[k] + len(substringa) < len(stringa) and (stringa[lpos[k] + len(substringa)].isdigit() == False and stringa[lpos[k] + len(substringa)].isalpha() == False and stringa[lpos[k] + len(substringa)] != '-'):
            stringa = stringa[:lpos[k]] + str(newStringa) + stringa[lpos[k] + len(substringa):]
            # lpos = [item + (lenNewStringa-len(substringa)) for item in lpos]
            lpos = [i.start() for i in re.finditer(substringa, stringa)]
            k = 0
            badPositions = 0
            # print('new lpos\t', lpos)
            # print('k update\t', k, '\t and badPositions\t', badPositions)
            # print('Update stringa\t', stringa, '\n')
        elif lpos[k] == 0 and lpos[k] + len(substringa) < len(stringa) and (stringa[lpos[k] + len(substringa)].isdigit() == False and \
            stringa[lpos[k] + len(substringa)].isalpha() == False and stringa[lpos[k] + len(substringa)] != '-'):
            stringa = stringa[:lpos[k]] + str(newStringa) + stringa[lpos[k] + len(substringa):]
            # lpos = [item + (lenNewStringa-len(substringa)) for item in lpos]
            lpos = [i.start() for i in re.finditer(substringa, stringa)]
            k = 0
            badPositions = 0
            # print('new lpos\t', lpos)
            # print('k update\t', k, '\t and badPositions\t', badPositions)
            # print('Update stringa\t', stringa, '\n')
        elif lpos[k] > 0 and lpos[k] + len(substringa) == len(stringa) and (stringa[lpos[k] - 1].isdigit() == False and \
            stringa[lpos[k] - 1].isalpha() == False and stringa[lpos[k] - 1] != '-'):
            stringa = stringa[:lpos[k]] + str(newStringa) + stringa[lpos[k] + len(substringa):]
            # lpos = [item + (lenNewStringa-len(substringa)) for item in lpos]
            lpos = [i.start() for i in re.finditer(substringa, stringa)]
            k = 0
            badPositions = 0
            # print('new lpos\t', lpos)
            # print('k update\t', k, '\t and badPositions\t', badPositions)
            # print('Update stringa\t', stringa, '\n')
        elif lpos[k] == 0 and lpos[k] + len(substringa) == len(stringa):
            stringa = str(newStringa)
            # lpos = [item + (lenNewStringa-len(substringa)) for item in lpos]
            lpos = [i.start() for i in re.finditer(substringa, stringa)]
            k = 0
            badPositions = 0
            # print('new lpos\t', lpos)
            # print('k update\t', k, '\t and badPositions\t', badPositions)
            # print('Update stringa\t', stringa, '\n')
        else:
            badPositions += 1
            k += 1
            # print('update badPositions\t', badPositions)
            # print('k update\t', k, '\t and badPositions\t', badPositions)
            # print('Update stringa\t', stringa, '\n')

    # print('Final stringa\t', stringa)
    return stringa

def replacePart(stringa, symbol):
    symbolPos = stringa.find(symbol)
    while symbolPos != -1:
        tmp = stringa[symbolPos:]
        spacePos = tmp.find(' ')
        if spacePos != -1:
            stringa = stringa[:symbolPos] + stringa[spacePos + symbolPos:]
        else:
            stringa = stringa[:symbolPos]
        symbolPos = stringa.find(symbol)
    return stringa

def getRandomObjectiveFunctions(nRxns,nFunctions):
    aOF = np.empty([nRxns,nFunctions])
    for i in range(aOF.shape[1]):
        allZeros = True
        while allZeros == True:
            tmpArray = np.random.random_sample((nRxns,))
            if np.all(tmpArray == 0) is True:
                allZeros = True
            else:
                allZeros = False
                aOF[:,i] = tmpArray
    return aOF

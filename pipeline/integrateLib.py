def computeScore(subrule, out, s, scoreAlreadyComputed):
    if ' or ' in subrule and ' and ' not in subrule:
        if out.empty is False and len(scoreAlreadyComputed) == 0:
            score = out[s].sum()
        if out.empty is False and len(scoreAlreadyComputed) != 0:
            score = out[s].sum() + sum(scoreAlreadyComputed)
        elif out.empty is True and len(scoreAlreadyComputed) == 0:
            score = float("NaN")
        elif out.empty is True and len(scoreAlreadyComputed) != 0:
            score = sum(scoreAlreadyComputed)
    elif ' or ' not in subrule and ' and ' in subrule:
        if out.empty is False and len(scoreAlreadyComputed) != 0:
            score = min(list(out[s]) + scoreAlreadyComputed)
        elif out.empty is False and len(scoreAlreadyComputed) == 0:
            score = min(list(out[s]))
        elif out.empty is True and len(scoreAlreadyComputed) == 0:
            score = float("NaN")
        elif out.empty is True and len(scoreAlreadyComputed) != 0:
            score = min(scoreAlreadyComputed)
    return score

def differenceKeepingDuplicates(v1, v2):
    outputDifference = []
    while v1:
        # Get first item (and remove).
        item = v1.pop(0)
        if item in v2:
            v2.remove(item)
        else:
            outputDifference.append(item)
    return outputDifference

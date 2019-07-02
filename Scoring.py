import sys


class OnTargetScore:

    def __init__(self):
        self.ScoringMatrix = []  # holds all the info of the Crisprscan features
        f = open('CASPERinfo')
        while True:
            line = f.readline()
            if line[:-1] == 'CRISPRSCAN_DATA':
                while True:
                    line = f.readline()
                    if line[0] == '-':
                        f.close()
                        break
                    myTup = line.split('\t')  # this is an input in the form nt1, nt2, pos, score
                    nt1 = myTup[0][0]
                    nt2 = myTup[0][1]
                    score = float(myTup[2][:-1])
                    iden = (nt1, nt2, int(myTup[1]), score)
                    self.ScoringMatrix.append(iden)
                break

    #  --- For every feature this def will check to see if it appears by checking the position in the sequence
    #  then confirm whether or not it is indeed the right nucleotide pattern --- #
    def returnScore(self, seq):
        totalScore = 0
        for iden in self.ScoringMatrix:
            pos = iden[2]
            if pos < len(seq):
                if iden[1] == 'x':
                    if seq[pos] == iden[0]:
                        totalScore += iden[3]
                else:
                    if seq[pos:pos+1] == iden[0]+iden[1]:
                        totalScore += iden[3]
        totalScore = 1-((1.29401-totalScore)/1.94947)
        totalScore = int((totalScore * 100) + 0.5)
        return totalScore


#dxy = OnTargetScore()
#print(dxy.returnScore("TTCAAAGTTCTGGGCAATAC"))
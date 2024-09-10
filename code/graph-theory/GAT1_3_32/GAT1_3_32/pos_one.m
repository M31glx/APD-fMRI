function Pos_One=pos_one(CorrMat)
CorrMat(CorrMat>0)=1;
Pos_One=CorrMat;
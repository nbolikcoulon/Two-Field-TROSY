#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import math
import numpy as np

import FitFunctions as FitF


def FitField(heights, fields):
    higherHeights = []
    lowerHeights = []
    middleHeights = []
    middleFields = []
    higherFields = []
    lowerFields = []
    for i in range(len(heights)):
        if heights[i] <= 0.47:
            higherHeights.append(heights[i])
            higherFields.append(fields[i])
        else:
            if heights[i] <= 0.6:
                middleHeights.append(heights[i])
                middleFields.append(fields[i])
            else:
                lowerHeights.append(heights[i])
                lowerFields.append(fields[i])
    
    HigherCoefs = np.polyfit(higherHeights, higherFields, 10)
    MiddleCoefs = np.polyfit(middleHeights, middleFields, 5)
    LowerCoefs = np.polyfit(lowerHeights, lowerFields, 10)
    
    return HigherCoefs, MiddleCoefs, LowerCoefs



#################################################### Fit of B0 ####################################################
def B0Fit(x, LowerCoefs, MiddleCoefs, HigherCoefs):
    Field = 0.0
    if x <= 0.47:
        for i in range(0, len(HigherCoefs)):
            Field += HigherCoefs[i] * math.pow(x, len(HigherCoefs) - i - 1)
    else:
        if x <= 0.6:
            for i in range(0, len(MiddleCoefs)):
                Field += MiddleCoefs[i] * math.pow(x, len(MiddleCoefs) - i - 1)
        else:
            for i in range(0, len(LowerCoefs)):
                Field += LowerCoefs[i] * math.pow(x, len(LowerCoefs) - i - 1)
    return Field




#################################################### Field List ####################################################
def FieldList(ShuttlingTime, Height, Increment, LowerCoefs, MiddleCoefs, HigherCoefs):
    Speed = Height/ShuttlingTime
    PositionList = np.arange(0, Height, Increment * Speed)
    PositionList = np.append(PositionList, Height)
    FieldList = [FitF.B0Fit(PositionList[i], LowerCoefs, MiddleCoefs, HigherCoefs) for i in range(len(PositionList))]
        

#    Acceleration = []
#    for i in range(len(self.ExperimentNumber)):
#        Acceleration.append(4.0*self.Height[i] / (self.SLF[i] * self.SLF[i]))
#    TimeListUp = [np.arange(0, self.SLF[i], Increment) for i in range(len(self.ExperimentNumber))]
#    TimeListDown = [np.arange(0, self.SHF[i], Increment) for i in range(len(self.ExperimentNumber))]
#    def getPosition(t, T, A, H):         #t for time, T for time to get to target point, A for acceleraiton, H target point
#        if t <= T/2.0:
#            h = 1.0 / 2.0 * A*t*t
#        else:
#            h = -1.0 / 2.0 * A * (t*t + T*T) + A*T*t + H
#        return h
#
#    PositionListUp = []
#    PositionListDown = []
#    for i in range(len(self.ExperimentNumber)):
#        IntermediateListUp = []
#        for Time in TimeListUp[i]:
#            IntermediateListUp.append(getPosition(Time, self.SLF[i], Acceleration[i], self.Height[i]))
#        IntermediateListDown = []
#        for Time in TimeListDown[i]:
#            IntermediateListDown.append(getPosition(Time, self.SHF[i], Acceleration[i], self.Height[i]))
#        PositionListUp.append(IntermediateListUp)
#        PositionListDown.append(IntermediateListDown)
#
#    FieldListUp = [[FitF.B0Fit(PositionListUp[i][j], self.LowerCoefs, self.MiddleCoefs, self.HigherCoefs) for j in range(len(PositionListUp[i]))] for i in range(len(PositionListUp))]
#    PositionListDown = [PositionListDown[i][::-1] for i in range(len(PositionListDown))]
#    FieldListDown = [[FitF.B0Fit(PositionListDown[i][j], self.LowerCoefs, self.MiddleCoefs, self.HigherCoefs) for j in range(len(PositionListDown[i]))] for i in range(len(PositionListDown))]
        
    return FieldList

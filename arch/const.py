class const:
    def __init__(self, finalBody, dvBudget, maxNumFlyby):
        self.dv = dvBudget
        self.flyby = maxNumFlyby
        self.finalBody = finalBody

    def check(self, finalBody, dv, numFlyby):
        constList = [finalBody != self.finalBody, dv < self.dv, numFlyby < self.flyby]
        return all(constList)
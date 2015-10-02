#!/usr/bin/python
from ConfigParser import *
import sys

class ConfigInformation:
    def __init__(self,Config):
        self._ListOfNecessarySections = ["CutScan","Signal","Background","Setup"]
        self._ListOfOptionalSection = ["Data"]
        self._SetupList = ["model","observable","theta","evaluationtype"]
        self.ListOfBkg  = []
        self.ListOfSig  = []
        self.ListOfCuts = []
        self.ModelDir = ""
        self.ThetaDir = ""
        self.Obs =""
        self.Eval ="" 
        self.Data = []
        self.fill_config(Config)

    def fill_config(self, Config):
        for option in self._SetupList:
            if not Config.has_option("Setup",option):
                print "Setup missing at least", option
                sys.exit(1)
        for necessarySection in self._ListOfNecessarySections:
            if necessarySection not in Config.sections():
                print necessarySection,"missing in the config File, can not procced"
                sys.exit(1)
            for option in Config.options(necessarySection):
                formated_option = [option,Config.get(necessarySection,option).split()]
                if necessarySection == "Background": self.ListOfBkg.append(formated_option)
                if necessarySection == "Signal": self.ListOfSig.append(formated_option)
                if necessarySection == "CutScan": self.ListOfCuts.append(formated_option)
                if necessarySection == "Setup":
                    if option == "model": self.ModelDir = Config.get(necessarySection,option)
                    if option == "observable": self.Obs = Config.get(necessarySection,option)
                    if option == "theta": self.ThetaDir = Config.get(necessarySection,option)
                    if option == "evaluationtype": self.Eval = Config.get(necessarySection,option)

        print "found necessary configuration"
        for optionalSection in self._ListOfOptionalSection:
            if optionalSection in Config.sections():
                for option in Config.options(optionalSection):
                    formated_option = [option]
                    formated_option.append(Config.get(necessarySection,option).split())
                    if optionalSection == "Data": Data.append(formated_option) 

    def PrintConfigInfo(self):
        print "Model Dir:          ", self.ModelDir
        print "Theta Dir:          ", self.ThetaDir
        print "Observable:         ", self.Obs
        print "Evaluatuion Typ:    ", self.Eval
        print "List of Backgrounds:", self.ListOfBkg
        print "List of Signals:    ", self.ListOfSig
        print "List of Cuts:       ", self.ListOfCuts
        if self.Data: print "Data Dir            ", self.Data


def ReadConfig(ConfigDir):
    Config = ConfigParser()
    Config.read(ConfigDir)
    print "starting to read configuration"
    return ConfigInformation(Config)

if __name__ == "__main__":
    TestInfo = ReadConfig("TestConfig.ini")
    TestInfo.PrintConfigInfo()
    


# Apr 23rd :-) 2019 trying to read movebank data using move

require(move)
#movebank user ors PW: eesH1eeF
MB.LoginObject=movebankLogin(username='ors',password='eesH1eeF')
#Gyps fulvus INPA Hatzofe Movebank ID	6071688
#HUJ MoveEcol Lab Israel: Griffon vulture Gyps fulvus Movebank ID	6638215
#timestamp_start, timestamp_end character or POSIXct=’yyyyMMddHHmmssSSS’
Dataset=getMovebankData(study=6071688, login=MB.LoginObject,
        includeExtraSensors=FALSE, deploymentAsIndividuals=FALSE,removeDuplicatedTimestamps=TRUE, ...)#animalName

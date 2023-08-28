# read data
dat <- haven::read_spss("//rtpnfil02/0272300.069_Wearables_CO_Innovation/Survey Data/Final files/MDOC Final Merged Database 7.28.23.sav")

# get unique variables
vnms <- names(dat)
vars <- unique(unlist(lapply(strsplit(vnms, "_"), \(x) x[1])))

# list daily log questions
dailies0 <- c(
"MindfulnessActivityCompletion",
"SleepQuality",
"DailySleepQuality",
"Work", 
#"WorkStart", 
#"WorkEnd", 
"WorkToday",
"WorkStress",
"WorkdayStress",
"Boredom", 
"Fatigue", 
"Overtime", 
"UndesirablePost", 
"ExtraWork", 
"AdminTasks", 
"InadequateStaffing", 
"DifficultyCoworker", 
"DifficultySupervisor", 
"InsufficientBreaks", 
"Threatened", 
"Assaulted", 
"WitnessedViolence", 
"NegativeNews", 
"Contraband", 
"DisciplinaryReport", 
"UsedForce", 
"OutsideWorkStress", 
"OutsideStress",
"MindfulnessActivity", 
"MindfulnessActivityCompletion",
"MindfulnessStress",
"MindfulnessStressRating")

# recreate daily dailies in order (44 max)
dailies <- lapply(dailies0, \(x){
  paste(x, 1:44, sep="_")
})

# scales


# non-daily data
nondailies <- names(dat)[! names(dat) %in% unlist(dailies)]


# make ID match what is in the wrap data
dat$id <- unlist(lapply(strsplit(dat$ParticipantID, "MDOCofficer"), \(x) as.numeric(x[2])))















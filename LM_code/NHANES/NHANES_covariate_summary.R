

data.path0 = './LM_code/NHANES/data/'
load(paste0(data.path0,'covariate_categorized_PA.Rdata'))

dim(Covariate)

### category to dummy

# Covariate = dummy_cols(Covariate,remove_selected_columns=TRUE,remove_most_frequent_dummy=TRUE)
# cols = colnames(Covariate)
# cols[16:20] = sapply(cols[16:20],function(x){gsub('.*_','',x)}) # marital(16-18), income(19-20)
# cols[21:22] = c('BMI.Normal','BMI.Obese') # BMI
# cols[23:26] = sapply(cols[23:26],function(x){gsub('.*_','',x)}) # race(23-25), gender(26)
# cols[27:31] = sapply(cols[27:31],function(x){gsub('_.*','',x)}) # diabetes(27), CHF(28), CHD(29), Cancer(30), Stroke(31)
# cols[32:35] = c('9-11th','Above.college','High.school','Less.9th') # education
# cols[36] = 'Mobility'
# cols[37:38] = sapply(cols[36:37],function(x){gsub('.*_','',x)}) # DrinkStatus
# cols[39:40] = c('Smoking','Smoking.Former') # smoking

Covariate = dummy_cols(Covariate,remove_selected_columns=TRUE,remove_most_frequent_dummy=TRUE)
cols = colnames(Covariate)
cols[19:23] = sapply(cols[19:23],function(x){gsub('.*_','',x)}) # marital(19-21), income(22-23)
cols[24:25] = c('BMI.Normal','BMI.Obese') # BMI
cols[26:29] = sapply(cols[26:29],function(x){gsub('.*_','',x)}) # race(26-28), gender(29)
cols[30:34] = sapply(cols[30:34],function(x){gsub('_.*','',x)}) # diabetes(30), CHF(31), CHD(32), Cancer(33), Stroke(34)
cols[35:38] = c('9-11th','Above.college','High.school','Less.9th') # education
cols[39] = 'Mobility'
cols[40:41] = sapply(cols[40:41],function(x){gsub('.*_','',x)}) # DrinkStatus
cols[42:43] = c('Smoking','Smoking.Former') # smoking

colnames(Covariate) = cols

dim(Covariate)

save(Covariate,file=paste0(data.path0,'Covariate_final.Rdata'))



# summary of variables

load(paste0(data.path0,'Covariate.Rdata'))

indices = Covariate$SEQN

load(paste0(data.path0,'Covariate_categorized.Rdata'))

Covariate = Covariate[Covariate$SEQN %in% indices,]


# summary of categorical variables
cols.category = c('Gender','Marital','Income','BMI_cat','Race','EducationAdult','DrinkStatus','SmokeCigs',
                  'Diabetes','CHF','CHD','Cancer','Stroke','MobilityProblem','Hypertensive')

n = nrow(Covariate)
print(n)
for (col in cols.category){
  print(col)
  print(rbind(table(Covariate[,col]),round(100*table(Covariate[,col])/n,2)))
}


# summary of continuous variables
cols.continuous = c('Age','BMI','Systolic.BP','Diastolic.BP','Height','Weight',
                    'Total.chol','HDL','LDL','Triglyceride','Glucose','Kcal','carbohydrate','protein','fat')
round(rbind(apply(Covariate[cols.continuous],2,function(x){c(mean(x),sd(x))})),2)


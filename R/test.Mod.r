### make this an s3 method for the anova generic

### test #!# Ty wants this to be KR, implement inside of pact.auto() as instructed
test.Mod <- function(mod1, mod2, data, KR=TRUE, KRmethod="ML", override=TRUE)
{
  p           <- NA
  test.Mod.err <- "test.Mod::Models are not nested"
  if(override)
  {
    #!# why is this override in here?? adding override override for PaCCT.Power.Nof1.r
    if(mod1$dims$Q==1 | mod1$dims$Q==1) KR <- FALSE
  }
  if("lme"%in%class(mod1) & "lme"%in%class(mod2))
  {
    mod.names   <- rbind(get.mod(mod1), get.mod(mod2))
    if(!KR)
    {
      try.lrt  <- try(mod.test <- anova(mod1, mod2), silent = TRUE)
      if(!"try-error"%in%is(try.lrt))
      {
        if('p-value'%in%names(mod.test))
        {
          mod.test <- data.frame(mod.names, mod.test)[,c(1:3, 5:12)]
          row.names(mod.test) <- NULL
          p <- mod.test$p.value[2]
        }
      }
      if("try-error"%in%is(try.lrt))
      {
        mod.test <- test.Mod.err
      }
    }
    if( KR)
    {
      mod1.r   <- lmeTOlmer(mod1, data, method=KRmethod)
      mod2.r   <- lmeTOlmer(mod2, data, method=KRmethod)
      try.kr <- try(mod.test <- KRmodcomp(mod1.r, mod2.r), silent = TRUE)
      if(!"try-error"%in%is(try.kr))
      {
        p <- mod.test$stats$p.value
      }
      if("try-error"%in%is(try.kr)) mod.test <- attr(try.kr, 'condition')
    }
  }
  # put warnings not captured by lmee in log
  #pact.err(warnings(), save_model, 'test.mod()')

  ### return tests
  return( list(test.type=ifelse(KR, "KR", "LRT"),
               models=mod.names,
               p=p))
}

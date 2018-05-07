# .Palytic - updates to the Palytic class

# this should only be applied to one participant at a time
Palytic$set("public", "getAR.order",
            function(w)
            {
              if(is.null(w)) stop('getAR.order() only works on 1 participant at a time')
              frmToChar(self$fixed)
              AR.order <- forecast::auto.arima(y)
              pq <- AR.order$arma[c(1,3)]


            },
            overwrite = TRUE)

Palytic$set("public", "getTime.Order",
            function(w, maxOrder=3)
            {
              if(is.null(w)) w <- 1:nrow(self$data)
              aics <- list()
              .v <- frmToChar(self$fixed)[1:2]
              .r <- frmToChar(self$random)
              for(i in 1:maxOrder)
              {
                .f <- formula( paste( .v[1], paste(.v[1], '^', i), collapse = '~' ) )
                aics[[i]] <- AIC( self$gamlss(fixed=.f, data=self$data, ) )
              }

            },
            overwrite = TRUE
)

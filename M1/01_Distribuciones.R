library(shiny)
library(ggplot2)

####################################################
# Paleta
####################################################

c1 <- "#813c60"
c2 <- "#cc6b96"
c3 <- "#be896d"
c4 <- "#ffe171"

####################################################
# UI
####################################################

ui <- fluidPage(
  
  titlePanel("Visualizador de Distribuciones"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      selectInput(
        "dist",
        "Distribución",
        choices = c(
          
          "Bernoulli",
          "Binomial",
          "Poisson",
          "Geométrica",
          "Binomial Negativa",
          
          "Uniforme",
          "Exponencial",
          "Gamma",
          "Beta",
          "Normal",
          "LogNormal",
          "Weibull",
          "Pareto",
          "t",
          "Ji-cuadrada",
          "F",
          "Cauchy",
          "Logística"
        )
      ),
      
      uiOutput("params")
      
    ),
    
    mainPanel(
      
      plotOutput("plot", height="650px")
      
    )
    
  )
  
)

####################################################
# SERVER
####################################################

server <- function(input, output){
  
  ####################################################
  # Parámetros dinámicos
  ####################################################
  
  output$params <- renderUI({
    
    switch(input$dist,
           
           "Bernoulli"=tagList(
             numericInput("p","p",0.5,0,1,0.01)
           ),
           
           "Binomial"=tagList(
             numericInput("n","n",20,1),
             numericInput("p","p",0.5,0,1,0.01)
           ),
           
           "Poisson"=tagList(
             numericInput("lambda","lambda",4,0.01)
           ),
           
           "Geométrica"=tagList(
             numericInput("p","p",0.30,0.001,1,0.01)
           ),
           
           "Binomial Negativa"=tagList(
             numericInput("size","size",5,1),
             numericInput("p","p",0.5,0.001,1,0.01)
           ),
           
           "Uniforme"=tagList(
             numericInput("min","mínimo",0),
             numericInput("max","máximo",10)
           ),
           
           "Exponencial"=tagList(
             numericInput("rate","rate",1)
           ),
           
           "Gamma"=tagList(
             numericInput("shape","shape",2),
             numericInput("rate","rate",1)
           ),
           
           "Beta"=tagList(
             numericInput("shape1","shape1",2),
             numericInput("shape2","shape2",5)
           ),
           
           "Normal"=tagList(
             numericInput("mean","media",0),
             numericInput("sd","desv",1,0.01)
           ),
           
           "LogNormal"=tagList(
             numericInput("meanlog","meanlog",0),
             numericInput("sdlog","sdlog",1)
           ),
           
           "Weibull"=tagList(
             numericInput("shape","shape",2),
             numericInput("scale","scale",1)
           ),
           
           "Pareto"=tagList(
             numericInput("shape","shape",3),
             numericInput("scale","xm",1)
           ),
           
           "t"=tagList(
             numericInput("df","gl",10,1)
           ),
           
           "Ji-cuadrada"=tagList(
             numericInput("df","gl",5)
           ),
           
           "F"=tagList(
             numericInput("df1","gl1",5),
             numericInput("df2","gl2",10)
           ),
           
           "Cauchy"=tagList(
             numericInput("location","location",0),
             numericInput("scale","scale",1)
           ),
           
           "Logística"=tagList(
             numericInput("location","location",0),
             numericInput("scale","scale",1)
           )
           
    )
    
  })
  
  ####################################################
  # Plot
  ####################################################
  
  output$plot <- renderPlot({
    
    ####################################################
    # DISCRETAS
    ####################################################
    
    if(input$dist=="Bernoulli"){
      
      x<-0:1
      y<-dbinom(x,1,input$p)
      
      ggplot(data.frame(x,y),aes(factor(x),y))+
        geom_col(fill=c2,width=.6)+
        theme_minimal(base_size=16)+
        labs(x="x",y="P(X=x)")
      
    }
    
    else if(input$dist=="Binomial"){
      
      x<-0:input$n
      y<-dbinom(x,input$n,input$p)
      
      ggplot(data.frame(x,y),aes(x,y))+
        geom_col(fill=c2)+
        theme_minimal(base_size=16)+
        labs(x = "x", y="P(X=x)")
      
    }
    
    else if(input$dist=="Poisson"){
      
      x<-0:max(20,qpois(.999,input$lambda))
      y<-dpois(x,input$lambda)
      
      ggplot(data.frame(x,y),aes(x,y))+
        geom_col(fill=c2)+
        theme_minimal(base_size=16)+
        labs(x = "x", y="P(X=x)")
      
    }
    
    else if(input$dist=="Geométrica"){
      
      x<-0:20
      y<-dgeom(x,input$p)
      
      ggplot(data.frame(x,y),aes(x,y))+
        geom_col(fill=c2)+
        theme_minimal(base_size=16)+
        labs(x = "x", y="P(X=x)")
      
    }
    
    else if(input$dist=="Binomial Negativa"){
      
      x<-0:30
      y<-dnbinom(x,input$size,input$p)
      
      ggplot(data.frame(x,y),aes(x,y))+
        geom_col(fill=c2)+
        theme_minimal(base_size=16)+
        labs(x = "x", y="P(X=x)")
      
    }
    
    ####################################################
    # CONTINUAS
    ####################################################
    
    else{
      
      if(input$dist=="Uniforme"){
        
        x<-seq(input$min,input$max,length=400)
        y<-dunif(x,input$min,input$max)
        
      }
      
      if(input$dist=="Exponencial"){
        
        x<-seq(0,qexp(.999,input$rate),length=400)
        y<-dexp(x,input$rate)
        
      }
      
      if(input$dist=="Gamma"){
        
        x<-seq(0,qgamma(.999,input$shape,input$rate),length=400)
        y<-dgamma(x,input$shape,input$rate)
        
      }
      
      if(input$dist=="Beta"){
        
        x<-seq(0,1,length=400)
        y<-dbeta(x,input$shape1,input$shape2)
        
      }
      
      if(input$dist=="Normal"){
        
        x<-seq(
          qnorm(.001,input$mean,input$sd),
          qnorm(.999,input$mean,input$sd),
          length=400)
        
        y<-dnorm(x,input$mean,input$sd)
        
      }
      
      if(input$dist=="LogNormal"){
        
        x<-seq(0,qlnorm(.999,input$meanlog,input$sdlog),length=400)
        y<-dlnorm(x,input$meanlog,input$sdlog)
        
      }
      
      if(input$dist=="Weibull"){
        
        x<-seq(0,qweibull(.999,input$shape,input$scale),length=400)
        y<-dweibull(x,input$shape,input$scale)
        
      }
      
      if(input$dist=="Pareto"){
        
        xm<-input$scale
        a<-input$shape
        
        x<-seq(xm,10*xm,length=400)
        
        y<-ifelse(x>=xm,a*xm^a/x^(a+1),0)
        
      }
      
      if(input$dist=="t"){
        
        x<-seq(qt(.001,input$df),qt(.999,input$df),length=400)
        y<-dt(x,input$df)
        
      }
      
      if(input$dist=="Ji-cuadrada"){
        
        x<-seq(0,qchisq(.999,input$df),length=400)
        y<-dchisq(x,input$df)
        
      }
      
      if(input$dist=="F"){
        
        x<-seq(0,qf(.999,input$df1,input$df2),length=400)
        y<-df(x,input$df1,input$df2)
        
      }
      
      if(input$dist=="Cauchy"){
        
        x<-seq(
          qcauchy(.001,input$location,input$scale),
          qcauchy(.999,input$location,input$scale),
          length=400)
        
        y<-dcauchy(x,input$location,input$scale)
        
      }
      
      if(input$dist=="Logística"){
        
        x<-seq(
          qlogis(.001,input$location,input$scale),
          qlogis(.999,input$location,input$scale),
          length=400)
        
        y<-dlogis(x,input$location,input$scale)
        
      }
      
      ggplot(data.frame(x,y),aes(x,y))+
        
        geom_area(fill=c1,alpha=.55)+
        
        geom_line(
          color=c1,
          linewidth=1.4
        )+
        
        theme_minimal(base_size=16)+
        
        theme(
          
          panel.grid.minor=element_blank(),
          
          plot.title=element_text(face="bold"),
          
          axis.title=element_text(color=c1)
          
        )+
        
        labs(
          
          x="x",
          y="Densidad"
          
        )
      
    }
    
  })
  
}

####################################################
# APP
####################################################

shinyApp(ui,server)
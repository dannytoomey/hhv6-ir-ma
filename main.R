lapply(c("googlesheets4","esc","meta","dplyr","metafor","grid",
    "brms","baggr","ggplot2","mclust"),require,character.only=TRUE)
setwd("./gosh_diagnostics")
files.sources = list.files()
sapply(files.sources, source)
setwd("../")

options(mc.cores = parallel::detectCores())
gs4_deauth()
options(error = function() traceback(3))

dataset <- read_sheet(
	"https://docs.google.com/spreadsheets/d/1qGOeQ-Jk_pvXYgzvlugtiVak8m9rKUmJDl8-dhoK1Qk/edit?gid=34534207#gid=34534207",
	sheet = "Included"
)
dataset <- dataset[dataset$`Included in analysis`=="Yes",]
save_path <- paste0("Data\ and\ Results/local_results/",gsub(":",";",substr(Sys.time(),1,19)),"")
if(!file.exists("Data\ and\ Results/local_results/")){
	dir.create("Data\ and\ Results/local_results/")
}
if(!file.exists(save_path)){
	dir.create(save_path)
}
if(!file.exists(paste0(save_path,"/plots"))){
	dir.create(paste0(save_path,"/plots"))
}
if(!file.exists(paste0(save_path,"/details"))){
	dir.create(paste0(save_path,"/details"))
}
neutrophil_text_save <- "/details/HHV-6 Detection and Neutrophil Engraftment.txt"
platelet_text_save <- "/details/HHV-6 Detection and Platelet Engraftment.txt"

if(!file.exists(paste0(save_path,neutrophil_text_save))){
    file.create(paste0(save_path,neutrophil_text_save))
}
if(!file.exists(paste0(save_path,platelet_text_save))){
    file.create(paste0(save_path,platelet_text_save))
}

or_analysis <- function(save_path,dataset,rma_pdf_width,bayes_pdf_width,
title,rma_pdf,rma_title,bayes_title,bayes_pdf,text_save,rma_pdf_height,
bayes_pdf_height,n_event_e,n_e,n_event_c,n_c,funnel_pdf,funnel_title){
    
    outlier_adjustment <- TRUE 

    cat(paste0("delayed engraftment ",Sys.time(),"\n"),
        file=paste0(save_path,text_save),append = TRUE)
    meta_bin <- metabin(
        event.e = unlist(dataset[c(n_event_e)],use.names=FALSE),
        n.e = unlist(dataset[c(n_e)],use.names=FALSE),
        event.c = unlist(dataset[c(n_event_c)],use.names=FALSE),
        n.c = unlist(dataset[c(n_c)],use.names=FALSE),
        studlab = dataset$`Study name`,
        sm = "OR",
        method = "MH",
        MH.exact = TRUE,
        fixed = FALSE,
        random = TRUE,
        method.tau = "PM",
        hakn = TRUE,
        title = title
    )

    if (outlier_adjustment){
        cat("\n-- RMA, OR, no outlier adjustment -- \n\n",
            file=paste0(save_path,text_save),append = TRUE)
        cat(paste0(capture.output(summary(meta_bin)),"\n"),
            file=paste0(save_path,text_save),append = TRUE)  

        meta_rma <- rma(
            yi = meta_bin$TE,
            sei = meta_bin$seTE,
            method = meta_bin$method.tau,
            test = "knha"
        )
        
        res_gosh <- gosh(meta_rma)
        res_gosh_diag <- gosh.diagnostics(res_gosh) 

        a <- res_gosh_diag$outlier.studies.km
        b <- res_gosh_diag$outlier.studies.db
        c <- res_gosh_diag$outlier.studies.gmm

        a[a > length(dataset$`Study name`)] <- length(dataset$`Study name`)
        b[b > length(dataset$`Study name`)] <- length(dataset$`Study name`)
        c[c > length(dataset$`Study name`)] <- length(dataset$`Study name`)

        if(length(a[a %in% b][a[a %in% b] %in% a[a %in% c]])>0){
            meta_bin <- update(meta_bin, exclude = c(a[a %in% b][a[a %in% b] %in% a[a %in% c]]))
        }

        cat("\n-- RMA, GOSH diagnostics -- \n\n",
            file=paste0(save_path,text_save),append = TRUE)
        cat(paste0(capture.output(print(res_gosh_diag)),"\n"),
            file=paste0(save_path,text_save),append = TRUE)
        cat("\n-- RMA, OR, outlier adjusted -- \n\n",
            file=paste0(save_path,text_save),append = TRUE)
        cat(paste0(capture.output(summary(meta_bin)),"\n"),
            file=paste0(save_path,text_save),append = TRUE)
    } else {
        cat("\n-- RMA, OR -- \n\n",
            file=paste0(save_path,text_save),append = TRUE)
        cat(paste0(capture.output(summary(meta_bin)),"\n"),
            file=paste0(save_path,text_save),append = TRUE)  
    }

    pdf(file = paste0(save_path,rma_pdf), width=rma_pdf_width, height=rma_pdf_height)
    metafor::forest(
        meta_bin,
        sortvar = TE,
        print.tau2 = TRUE,
        leftlabs = c("Study", "N delayed","Total","N delayed","Total"),
        lab.e = "HHV-6 +",
        lab.c = "HHV-6 -"
    )
    grid.text(rma_title, x=0.5,y=0.95, gp=gpar(fontsize=16))

    if(dim(dataset)[1] >= 10){
        # use peters bias for binary data
        cat("\n-- Publication Bias -- \n\n",file=paste0(save_path,text_save),append = TRUE)
        cat(paste0(capture.output(metabias(meta_bin, method.bias = "peters")),"\n"),
            file=paste0(save_path,text_save),append = TRUE)
        pdf(file = paste0(save_path,funnel_pdf), width = 10, height = 10)
        meta::funnel(meta_bin,studlab=TRUE)
        grid.text(funnel_title, x=0.5,y=0.95, gp=gpar(fontsize=16))
        dev.off()
    }
    
    if (outlier_adjustment){
        prep_ma_df <- data.frame(group = meta_bin$studlab, tau = meta_bin$TE, 
            se = meta_bin$seTE, exclude=meta_bin$exclude)
        prep_ma_df <- prep_ma_df[prep_ma_df$exclude==FALSE,]
    } else {
        prep_ma_df <- data.frame(group = meta_bin$studlab, tau = meta_bin$TE, 
            se = meta_bin$seTE)
    }
    baggr_bin <- baggr(prep_ma_df, model = "rubin", pooling = "partial",
        iter=25000, chains=10)
    if (outlier_adjustment){
        cat("\n-- Bayes, OR, outlier adjusted -- \n",
            file=paste0(save_path,text_save),append = TRUE)
    } else {
        cat("\n-- Bayes, OR -- \n",file=paste0(save_path,text_save),append = TRUE)
    }

    if (max(rstan:::summary_sim(baggr_bin$fit@sim)$rhat) > 1.05) {
        cat("Operation terminated due to lack of chain convergence. \n",
            file=paste0(save_path,text_save),append = TRUE)
    } else {
        # save model fit and results if chains converged
        cat("\n- Model fit - \n",file=paste0(save_path,text_save),append = TRUE)
        cat(paste0(capture.output(print(baggr_bin$fit)),"\n"),
            file=paste0(save_path,text_save),append = TRUE)
        cat("\n- Model results - \n",file=paste0(save_path,text_save),append = TRUE)
        cat(paste0(capture.output(print(baggr_bin)),"\n"),
            file=paste0(save_path,text_save),append = TRUE)
        # save baggr plot
        p <- plot(baggr_bin,hyper=TRUE,vline=FALSE,style="areas") + bayes_title
        pdf(file = paste0(save_path,bayes_pdf), width=bayes_pdf_width, 
            height=bayes_pdf_height, onefile=FALSE)
        print(p)
        dev.off()
    }    
}

neutrophil_dataset <- dataset[unlist(lapply(dataset$`N HHV-6 pos pts with delayed or failed neutrophil engraftment`,is.numeric)) ,] 
rma_pdf_width <- 10
bayes_pdf_width <- 8
neutrophil_title <- "HHV-6 and Delayed Neutrophil Engraftment"
neutrophil_rma_pdf <- "/plots/RMA OR HHV-6 Delayed Neutrophil Engraftment.pdf"
neutrophil_funnel_pdf <- "/plots/Funnel OR HHV-6 Delayed Neutrophil Engraftment.pdf"
neutrophil_funnel_title <- "Publication Bias for HHV-6 and Delayed Neutrophil Engraftment"
neutrophil_rma_title <- "Risk of Delayed Neutrophil Engraftment with HHV-6 positivity"
neutrophil_bayes_title <- ggtitle("Bayesian Aggregations: Delayed Neutrophil Engraftment and HHV-6 positivity",
    "Posterior distributions with 95% intervals")
neutrophil_bayes_pdf <- "/plots/Bayes OR HHV-6 Delayed Neutrophil Engraftment.pdf"
neutrophil_rma_pdf_height <- log(dim(neutrophil_dataset)[1],9)*5.25
neutrophil_bayes_pdf_height <- log(dim(neutrophil_dataset)[1],9)*5.25
neutrophil_n_event_e <- "N HHV-6 pos pts with delayed or failed neutrophil engraftment"
neutrophil_n_e <- "N patients with HHV-6"
neutrophil_n_event_c <- "N HHV-6 neg pts with delayed or failed neutrophil engraftment"
neutrophil_n_c <- "N patients without HHV-6"

or_analysis(save_path,neutrophil_dataset,rma_pdf_width,bayes_pdf_width,neutrophil_title,
neutrophil_rma_pdf,neutrophil_rma_title,neutrophil_bayes_title,neutrophil_bayes_pdf,
neutrophil_text_save,neutrophil_rma_pdf_height,neutrophil_bayes_pdf_height,
neutrophil_n_event_e,neutrophil_n_e,neutrophil_n_event_c,neutrophil_n_c,
neutrophil_funnel_pdf,neutrophil_funnel_title)

platelet_dataset <- dataset[!is.na(dataset$`N HHV-6 pos pts with delayed or failed platelet engraftment`) ,] 
platelet_title <- "HHV-6 and Delayed Platelet Engraftment"
platelet_rma_pdf <- "/plots/RMA OR HHV-6 Delayed Platelet Engraftment.pdf"
platelet_funnel_pdf <- "/plots/Funnel OR HHV-6 Delayed Platelet Engraftment.pdf"
platelet_funnel_title <- "Publication Bias for HHV-6 and Delayed Platelet Engraftment"
platelet_rma_title <- "Risk of Delayed Platelet Engraftment with HHV-6 positivity"
platelet_bayes_title <- ggtitle("Bayesian Aggregations: Delayed Platelet Engraftment and HHV-6 positivity",
    "Posterior distributions with 95% intervals")
platelet_bayes_pdf <- "/plots/Bayes OR HHV-6 Delayed Platelet Engraftment.pdf"
platelet_rma_pdf_height <- log(dim(platelet_dataset)[1],9)*5.25
platelet_bayes_pdf_height <- log(dim(platelet_dataset)[1],9)*5.25
platelet_n_event_e <- "N HHV-6 pos pts with delayed or failed platelet engraftment"
platelet_n_e <- "N patients with HHV-6"
platelet_n_event_c <- "N HHV-6 neg pts with delayed or failed platelet engraftment"
platelet_n_c <- "N patients without HHV-6"

or_analysis(save_path,platelet_dataset,rma_pdf_width,bayes_pdf_width,platelet_title,
platelet_rma_pdf,platelet_rma_title,platelet_bayes_title,platelet_bayes_pdf,
platelet_text_save,platelet_rma_pdf_height,platelet_bayes_pdf_height,
platelet_n_event_e,platelet_n_e,platelet_n_event_c,platelet_n_c,
platelet_funnel_pdf,platelet_funnel_title)

if (file.exists("./Rplots.pdf")) {
  file.remove("./Rplots.pdf")
}

# https://training.cochrane.org/sites/training.cochrane.org/files/public/uploads/resources/downloadable_resources/meta-analysis%20of%20time%20to%20event%20data%20Cochrane%20webinar%20July%202018.pdf
# https://www.pharmasug.org/proceedings/2023/DV/PharmaSUG-2023-DV-287.pdf
# https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html
# https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781444311723.oth2

lapply(c("googlesheets4","esc","meta","dplyr","metafor","grid","brms","baggr","ggplot2","dmetar"),require,character.only=TRUE)
library("metamedian",lib.loc="./metamedian_dev/")

options(mc.cores = parallel::detectCores())
gs4_deauth()
options(error = function() traceback(3))

dataset <- read_sheet(
	"https://docs.google.com/spreadsheets/d/1qLUi0n5yuLv9DyjAVRSfr2NmtdO2N6wvFutUYKfYnZ0/edit#gid=1287100954",
	sheet = "Included"
)
dataset <- dataset[dataset$`Included in analysis`=="Yes",]
save_path <- paste0("local_results/",gsub(":",";",substr(Sys.time(),1,19)),"")
if(!file.exists("local_results/")){
	dir.create("local_results/")
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

or_analysis <- function(save_path,dataset,rma_pdf_width,bayes_pdf_width,title,rma_pdf,rma_title,bayes_title,bayes_pdf,text_save,rma_pdf_height,bayes_pdf_height,n_event_e,n_e,n_event_c,n_c){
    cat(paste0("delayed engraftment ",Sys.time(),"\n"),file=paste0(save_path,text_save),append = TRUE)
    meta_bin <- metabin(
        event.e = unlist(dataset[c(n_event_e)],use.names=FALSE),
        n.e = unlist(dataset[c(n_e)],use.names=FALSE),
        event.c = unlist(dataset[c(n_event_c)],use.names=FALSE),
        n.c = unlist(dataset[c(n_c)],use.names=FALSE),
        studlab = dataset$`Study name`,
        sm = "RR",
        method = "MH",
        MH.exact = TRUE,
        fixed = FALSE,
        random = TRUE,
        method.tau = "PM",
        hakn = TRUE,
        title = title
    )
    cat("\n-- RMA, RR, no outlier adjustment -- \n\n",file=paste0(save_path,text_save),append = TRUE)
    cat(paste0(capture.output(summary(meta_bin)),"\n"),file=paste0(save_path,text_save),append = TRUE)
    meta_rma <- rma(
        yi = meta_bin$TE,
        sei = meta_bin$seTE,
        method = meta_bin$method.tau,
        test = "knha"
    )
    res_gosh <- gosh(meta_rma)
    res_gosh_diag <- gosh.diagnostics(res_gosh)
    a<-res_gosh_diag$outlier.studies.km
    b<-res_gosh_diag$outlier.studies.db
    c<-res_gosh_diag$outlier.studies.gmm
    if(length(a[a %in% b][a[a %in% b] %in% a[a %in% c]])>0){
        meta_bin <- update(meta_bin, exclude = c(a[a %in% b][a[a %in% b] %in% a[a %in% c]]))
    }
    pdf(file = paste0(save_path,rma_pdf), width=rma_pdf_width, height=rma_pdf_height)
    meta::forest.meta(
        meta_bin,
        sortvar = TE,
        print.tau2 = TRUE,
        leftlabs = c("Study", "N delayed","Total","N delayed","Total"),
        lab.e = "HHV-6 +",
        lab.c = "HHV-6 -"
    )
    grid.text(rma_title, x=0.5,y=0.95, gp=gpar(fontsize=16))

    prep_ma_df <- data.frame(group = meta_bin$studlab, tau = meta_bin$TE, se = meta_bin$seTE, exclude=meta_bin$exclude)
    prep_ma_df <- prep_ma_df[prep_ma_df$exclude==FALSE,]
    baggr_bin <- baggr(prep_ma_df, model = "rubin", pooling = "partial",iter=20000,chains=10)
    p <- plot(baggr_bin,hyper=TRUE,vline=FALSE,style="areas") + bayes_title
    pdf(file = paste0(save_path,bayes_pdf), width=bayes_pdf_width, height=bayes_pdf_height, onefile=FALSE)
    print(p)
    dev.off()

    cat("\n-- RMA, GOSH diagnostics -- \n\n",file=paste0(save_path,text_save),append = TRUE)
    cat(paste0(capture.output(print(res_gosh_diag)),"\n"),file=paste0(save_path,text_save),append = TRUE)
    cat("\n-- RMA, RR, outlier adjusted -- \n\n",file=paste0(save_path,text_save),append = TRUE)
    cat(paste0(capture.output(summary(meta_bin)),"\n"),file=paste0(save_path,text_save),append = TRUE)
    cat("\n-- Bayes, RR -- \n\n",file=paste0(save_path,text_save),append = TRUE)
    cat(paste0(capture.output(print(baggr_bin)),"\n"),file=paste0(save_path,text_save),append = TRUE)

}

analysis <- c("RR")

for(a in analysis){
    if(a == "RR"){
        neutrophil_dataset <- dataset[unlist(lapply(dataset$`N HHV-6 pos pts with delayed or failed neutrophil engraftment`,is.numeric)) & dataset$`N patients with HHV-6`>10,]
        rma_pdf_width <- 10
        bayes_pdf_width <- 8
        neutrophil_title <- "HHV-6 and Delayed Neutrophil Engraftment"
        neutrophil_rma_pdf <- "/plots/RMA RR HHV-6 Delayed Neutrophil Engraftment.pdf"
        neutrophil_rma_title <- "Risk of Delayed Neutrophil Engraftment with HHV-6 positivity"
        neutrophil_bayes_title <- ggtitle("Bayesian Aggregations: Delayed Neutrophil Engraftment and HHV-6 positivity","Posterior distributions with 95% intervals")
        neutrophil_bayes_pdf <- "/plots/Bayes RR HHV-6 Delayed Neutrophil Engraftment.pdf"
        neutrophil_rma_pdf_height <- log(dim(neutrophil_dataset)[1],9)*5.25
        neutrophil_bayes_pdf_height <- log(dim(neutrophil_dataset)[1],9)*5.25
        neutrophil_n_event_e <- "N HHV-6 pos pts with delayed or failed neutrophil engraftment"
        neutrophil_n_e <- "N patients with HHV-6"
        neutrophil_n_event_c <- "N HHV-6 neg pts with delayed or failed neutrophil engraftment"
        neutrophil_n_c <- "N patients without HHV-6"

        or_analysis(save_path,neutrophil_dataset,rma_pdf_width,bayes_pdf_width,neutrophil_title,
        neutrophil_rma_pdf,neutrophil_rma_title,neutrophil_bayes_title,neutrophil_bayes_pdf,
        neutrophil_text_save,neutrophil_rma_pdf_height,neutrophil_bayes_pdf_height,
        neutrophil_n_event_e,neutrophil_n_e,neutrophil_n_event_c,neutrophil_n_c)

        platelet_dataset <- dataset[!is.na(dataset$`N HHV-6 pos pts with delayed or failed platelet engraftment`) & dataset$`N patients with HHV-6`>10,]
        platelet_title <- "HHV-6 and Delayed Platelet Engraftment"
        platelet_rma_pdf <- "/plots/RMA RR HHV-6 Delayed Platelet Engraftment.pdf"
        platelet_rma_title <- "Risk of Delayed Platelet Engraftment with HHV-6 positivity"
        platelet_bayes_title <- ggtitle("Bayesian Aggregations: Delayed Platelet Engraftment and HHV-6 positivity","Posterior distributions with 95% intervals")
        platelet_bayes_pdf <- "/plots/Bayes RR HHV-6 Delayed Platelet Engraftment.pdf"
        platelet_rma_pdf_height <- log(dim(platelet_dataset)[1],9)*5.25
        platelet_bayes_pdf_height <- log(dim(platelet_dataset)[1],9)*5.25
        platelet_n_event_e <- "N HHV-6 pos pts with delayed or failed platelet engraftment"
        platelet_n_e <- "N patients with HHV-6"
        platelet_n_event_c <- "N HHV-6 neg pts with delayed or failed platelet engraftment"
        platelet_n_c <- "N patients without HHV-6"

        or_analysis(save_path,platelet_dataset,rma_pdf_width,bayes_pdf_width,platelet_title,
        platelet_rma_pdf,platelet_rma_title,platelet_bayes_title,platelet_bayes_pdf,
        platelet_text_save,platelet_rma_pdf_height,platelet_bayes_pdf_height,
        platelet_n_event_e,platelet_n_e,platelet_n_event_c,platelet_n_c)
    }
    if(a == "median"){
        neutro_meta_med_range <- metamedian(
            data = dataset[!is.na(dataset$`HHV-6 pos min days to neutrophil engraftment`),],
           	median_method = "qe",
           	group_1_n="N patients with HHV-6",
           	group_1_min="HHV-6 pos min days to neutrophil engraftment",
           	group_1_med="HHV-6 pos median days to neutrophil engraftment",
           	group_1_max="HHV-6 pos max days to neutrophil engraftment",
           	group_2_n="N patients without HHV-6",
           	group_2_min="HHV-6 neg min days to neutrophil engraftment",
           	group_2_med="HHV-6 neg median days to neutrophil engraftment",
           	group_2_max="HHV-6 neg max days to neutrophil engraftment"
        )
        platelet_meta_med_range <- metamedian(
            data = dataset[!is.na(dataset$`HHV-6 neg min days to platelet engraftment`),],
           	median_method = "qe",
           	group_1_n="N patients with HHV-6",
           	group_1_min="HHV-6 pos min days to platelet engraftment",
           	group_1_med="HHV-6 pos median days to platelet engraftment",
           	group_1_max="HHV-6 pos max days to platelet engraftment",
           	group_2_n="N patients without HHV-6",
           	group_2_min="HHV-6 neg min days to platelet engraftment",
           	group_2_med="HHV-6 neg median days to platelet engraftment",
           	group_2_max="HHV-6 neg max days to platelet engraftment"
        )
        cat("-- RMA, median and range -- \n",file=paste0(save_path,neutrophil_text_save),append = TRUE)
        cat(paste0(capture.output(print(neutro_meta_med_range)),"\n"),file=paste0(save_path,neutrophil_text_save),append = TRUE)

        cat("-- RMA, median and range -- \n",file=paste0(save_path,platelet_text_save),append = TRUE)
        cat(paste0(capture.output(print(platelet_meta_med_range)),"\n"),file=paste0(save_path,platelet_text_save),append = TRUE)
    }
    if(a == "HR"){
        neutrophil_data = dataset[!is.na(dataset$`Combined HR/MHR neutrophil engraftment`),]
        platelet_data = dataset[!is.na(dataset$`Combined HR/MHR platelet engraftment`),]

        neutrophil_meta<-metagen(
            HR = log(neutrophil_data$`Combined HR/MHR neutrophil engraftment`),
            lower = log(neutrophil_data$`Combined HR/MHR LB neutrophil engraftment`),
            upper = log(neutrophil_data$`Combined HR/MHR UB neutrophil engraftment`),
            sm = "HR", fixed=F, random=T,
            studlab = neutrophil_data$`Study name`,
            method.tau = "DL",
            method.random.ci = "classic",
        )
        platelet_meta<-metagen(
            HR = log(platelet_data$`Combined HR/MHR platelet engraftment`),
            lower = log(platelet_data$`Combined HR/MHR LB platelet engraftment`),
            upper = log(platelet_data$`Combined HR/MHR UB platelet engraftment`),
            sm = "HR", fixed=F, random=T,
            studlab = platelet_data$`Study name`,
            method.tau = "DL",
            method.random.ci = "classic",
        )

        rma_pdf_width <- 10
        neutrophil_rma_pdf <- "/plots/RMA HR HHV-6 Neutrophil Engraftment.pdf"
        neutrophil_rma_pdf_height <- log(dim(neutrophil_data)[1],9)*5.25
        neutrophil_rma_title <- "HR of Neutrophil Engraftment with HHV-6 positivity"
        pdf(file = paste0(save_path,neutrophil_rma_pdf), width=rma_pdf_width, height=neutrophil_rma_pdf_height)
        meta::forest.meta(
            neutrophil_meta,
            sortvar = TE,
            print.tau2 = TRUE,
        )
        grid.text(neutrophil_rma_title, x=0.5,y=0.95, gp=gpar(fontsize=16))

        platelet_rma_pdf <- "/plots/RMA HR HHV-6 Platelet Engraftment.pdf"
        platelet_rma_pdf_height <- log(dim(platelet_data)[1],9)*5.25
        platelet_rma_title <- "HR of Platelet Engraftment with HHV-6 positivity"
        pdf(file = paste0(save_path,platelet_rma_pdf), width=rma_pdf_width, height=platelet_rma_pdf_height)
        meta::forest.meta(
            platelet_meta,
            sortvar = TE,
            print.tau2 = TRUE,
        )
        grid.text(platelet_rma_title, x=0.5,y=0.95, gp=gpar(fontsize=16))

        neutrophil_baggr_df <- data.frame(group = neutrophil_meta$studlab, tau = neutrophil_meta$TE, se = neutrophil_meta$seTE)
        platelet_baggr_df <- data.frame(group = platelet_meta$studlab, tau = platelet_meta$TE, se = platelet_meta$seTE)

        bayes_pdf_width <- 8
        neutrophil_bayes_title <- ggtitle("Hazard Ratio of Time to Neutrophil Engraftment and HHV-6 positivity","Posterior distributions with 95% intervals")
        neutrophil_bayes_pdf <- "/plots/Bayes HR HHV-6 Neutrophil Engraftment.pdf"
        neutrophil_bayes_pdf_height <- log(dim(neutrophil_dataset)[1],9)*5.25
        neutrophil_baggr_bin <- baggr(neutrophil_baggr_df, model = "rubin", pooling = "partial",iter=20000,chains=10)
        p <- plot(neutrophil_baggr_bin,hyper=TRUE,vline=FALSE,style="areas") + neutrophil_bayes_title
        pdf(file = paste0(save_path,neutrophil_bayes_pdf), width=bayes_pdf_width, height=neutrophil_bayes_pdf_height, onefile=FALSE)
        print(p)
        dev.off()

        platelet_bayes_title <- ggtitle("Hazard Ratio of Time to Platelet Engraftment and HHV-6 positivity","Posterior distributions with 95% intervals")
        platelet_bayes_pdf <- "/plots/Bayes HR HHV-6 Platelet Engraftment.pdf"
        platelet_bayes_pdf_height <- log(dim(neutrophil_dataset)[1],9)*5.25
        platelet_baggr_bin <- baggr(platelet_baggr_df, model = "rubin", pooling = "partial",iter=20000,chains=10)
        p <- plot(platelet_baggr_bin,hyper=TRUE,vline=FALSE,style="areas") + platelet_bayes_title
        pdf(file = paste0(save_path,platelet_bayes_pdf), width=bayes_pdf_width, height=platelet_bayes_pdf_height, onefile=FALSE)
        print(p)
        dev.off()

        cat("-- RMA, HR and MHR -- \n\n",file=paste0(save_path,neutrophil_text_save),append = TRUE)
        cat(paste0(capture.output(print(summary(neutrophil_meta))),"\n"),file=paste0(save_path,neutrophil_text_save),append = TRUE)
        cat("\n-- Bayes, HR and MHR -- \n\n",file=paste0(save_path,neutrophil_text_save),append = TRUE)
        cat(paste0(capture.output(print(neutrophil_baggr_bin)),"\n"),file=paste0(save_path,neutrophil_text_save),append = TRUE)

        cat("-- RMA, HR and MHR -- \n\n",file=paste0(save_path,platelet_text_save),append = TRUE)
        cat(paste0(capture.output(print(summary(platelet_meta))),"\n"),file=paste0(save_path,platelet_text_save),append = TRUE)
        cat("\n-- Bayes, HR and MHR -- \n\n",file=paste0(save_path,platelet_text_save),append = TRUE)
        cat(paste0(capture.output(print(platelet_baggr_bin)),"\n"),file=paste0(save_path,platelet_text_save),append = TRUE)
    }
}

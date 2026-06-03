
rm(list=ls())

library(ggplot2)
library(stringr)
#install.packages("patchwork")
#library(patchwork)
#install.packages("staplr")
library(staplr)
library(readxl)

setwd("Z:\\EMA_causal_16May_node2\\CI.RCT.Sim\\results_24Mai")
files<-list.files(recursive=TRUE)
#files

#setwd("Z:\\EMA_causal_16May_node2\\CI.RCT.Sim\\results4")

#load("ALL_logHR_20May_1cscen_small_n_H1_sc_1to8_results_onco_nsim100_ims-node22026-05-20_1446.Rdata")
#results$cens.est.coverage
#load("ALL_logHR_19May_Mscen_small_n_H1_sc_1to8_results_onco_nsim1000_ims-node132026-05-19_2227.RData")

#setwd("C:\\EMA_Causal\\CI.RCT.Sim\\results3")
scenario_list <- as.data.frame(read_xlsx("oncology_scenario_list.xlsx"))
#correct the name of progression - fast to progression - slow
scenario_list$scenario_name<-gsub("progression - fast","progression - slow",scenario_list$scenario_name)


get_pos<-function(x,pos=1) {
		spl<-strsplit(x,"[.]")
		out<-rep(NA,length(spl))
		for(i in 1:length(spl)) out[i]<-spl[[i]][pos]
		out
}


#names(results)
make_tab<-function(file_sel) {
	read_files<-files[file_sel]

	stat_sel<-c("est.bias","est.sd_bias","rejection_0.025","est.coverage","est.width","est.mse")
	stat_name<-c("Bias",   "SE_bias",    "Rejection_rate", "CI_coverage", "CI_width",  "RMSE")

	#stat_sel<-"est.bias"
	all_tab<-vector(length=length(stat_sel)+1,mode="list") #+1 for descriptive summary statistics
	#i<-1
	for(i in 1:length(read_files)) {
		print(i)
		flush.console()
		load(read_files[i])
		#results structure
		##scenario, iteration, methode, parameter
		results<-as.data.frame(results)
		#j<-1
		for(j in 1:length(stat_sel)) {
			ind<-grepl(stat_sel[j],names(results))
			tab0<-results[,ind]
			#
			#results[1,1,"
			#
			names(tab0)<-get_pos(names(tab0),1)
			tab1<-data.frame(nsim=results$REPLICATIONS,statistic=stat_sel[j],scenario=results$scen_set,tab0)
			if(stat_sel[j]=="est.sd_bias") {
				tab1[,-(1:3)]<-tab1[,-(1:3)]/sqrt(tab1$nsim)
				tab1$statistic<-"est.SE_bias"
			}
			if(stat_sel[j]=="est.mse") {
				tab1[,-(1:3)]<-sqrt(tab1[,-(1:3)])
				tab1$statistic<-"RMSE"
			}

			all_tab[[j]]<-rbind(all_tab[[j]],tab1)
		}
		#descriptives

		r0<-function(x) round(x,0)
		r<-function(x,a=1) formatC(x,digits=a,format="f")
		tab_descr<-with(results,data.frame(
			#Scenario=scenario_list$scenario_name[match(scen_set,scenario_list$Scen_ID)],
			scenario=scen_set,
			n_trt=r(describe.n_trt),
			n_ctr=r(describe.n_ctrl),
			ev_trt=r(describe.ev_trt),
			ev_ctr=r(describe.ev_ctrl),
			n_switch=r(describe.n_switch),
			n_rcens=r(describe.sd_n_random_cens),
			ev_achieved=r(describe.sufficient_events,4),  #*REPLICATIONS
			max_FU_years=r(describe.max_followup),
			trueHR=r(exp(true_eff),a=2)
		))
		j<-j+1
		#all_tab[[j]]<-tab_descr
		all_tab[[j]]<-rbind(all_tab[[j]],tab_descr)

		#all_tab[[i]]<-tab
	}
	all_tab

	i<-1
	for(i in 1:length(all_tab)) {
		all_tab[[i]]<-all_tab[[i]][order(all_tab[[i]]$scenario),]
		all_tab[[i]]<-merge(all_tab[[i]],scenario_list[,1:3],by.x="scenario",by.y="Scen_ID",all.x=TRUE,all.y=FALSE)
	}
	names(all_tab)<-c(stat_name,"Descriptives")
	all_tab
}



name_small<-"Sim_25Mai_korr_10kscen" ###"results_24Mai/Sim_24May_1c"
name_large<-"Sim_25Mai_korr_10kscen" ###"results4/ALL_logHR_20May_10k"
name2<-"nsim100"

tab_list<-list(
	small_n_H1=make_tab(file_sel=grepl(name_small,files)&grepl(name2,files)&grepl("small_n",files)&grepl("H1",files) ),
	small_n_H0=make_tab(grepl(name_small,files)&grepl(name2,files)&grepl("small_n",files)&grepl("H0",files) ),

	large_n_H1=make_tab(grepl(name_large,files)&grepl(name2,files)&grepl("large_n",files)&grepl("H1",files) ),
	large_n_H0=make_tab(grepl(name_large,files)&grepl(name2,files)&grepl("large_n",files)&grepl("H0",files) ),

	extra_cens_H1=make_tab(grepl(name_small,files)&grepl(name2,files)&grepl("IPCW_extra",files)&grepl("H1",files) ),
	extra_cens_H0=make_tab(grepl(name_small,files)&grepl(name2,files)&grepl("IPCW_extra",files)&grepl("H0",files) )
)


library(writexl)


if(!file.exists("Grafiken_27Mai"))  dir.create("Grafiken_27Mai")

#setwd("/home/robin/EMA_causal_16May_node2/CI.RCT.Sim/results_24Mai/Grafiken_27Mai")
setwd("Z:\\EMA_causal_16May_node2\\CI.RCT.Sim\\results_24Mai\\Grafiken_27Mai")

save(tab_list,file="tab_list.RData")


#save descriptive tabs to latex
i<-1
for(i in 1:length(tab_list)) {

	write_xlsx(tab_list[[i]],paste(names(tab_list)[i],"_results.xlsx",sep=""))
}

filename<-paste("Descriptives.txt",sep="")
i<-1
for(i in 1:length(tab_list)) {
	savtab<-tab_list[[i]]$Descriptives

	#savtab$Scen_ID<-savtab$scenario
	savtab$scenario<-savtab$scenario_name
	savtab$scenario_name<-NULL
	savtab$block_nam<-NULL
	names(savtab)<-gsub("n_trt","trt",names(savtab))
	names(savtab)<-gsub("n_ctr","ctr",names(savtab))
	names(savtab)<-gsub("n_switch","switch",names(savtab))
	names(savtab)<-gsub("n_rcens","rcens",names(savtab))
	savtab$ev_achieved<-as.character(round((1-as.numeric(savtab$ev_achieved))*tab_list[[i]]$Bias$nsim[1])) #round to correct for numerical issues within 1-small number
	names(savtab)<-gsub("ev_achieved","ev<aim",names(savtab))
	names(savtab)<-gsub("max_FU_years","FU",names(savtab))
	names(savtab)<-gsub("trueHR","HR",names(savtab))

	names(savtab)<-gsub("_"," ",names(savtab))

	blocknam<-gsub("_"," ",names(tab_list)[i])
	if(blocknam=="small n H0") blocknam<-"scenarios with small sample size under the null hypothesis."
	if(blocknam=="small n H1") blocknam<-"scenarios with small sample size under the alternative hypothesis."

	if(blocknam=="large n H0") blocknam<-"scenarios with large sample size under the null hypothesis."
	if(blocknam=="large n H1") blocknam<-"scenarios with large sample size under the alternative hypothesis."

	if(blocknam=="extra cens H0") blocknam<-"scenarios with different random censoring patterns under the null hypothesis."
	if(blocknam=="extra cens H1") blocknam<-"scenarios with different random censoring patterns under the alternatve hypothesis."



	#savtab<-data.frame(Block=blocknam,savtab)
	savtab
	#filename<-paste(names(tab_list)[i],"_Descriptives.txt",sep="")

	x<-savtab$scenario[4]
	fun<-function(x,w=40) {
		nc<-nchar(x)
		if(nc>w) {
			split<-ceiling(nc/2)
			x<-paste(substr(x,1,split),"MMM",substr(x,split+1,nc),sep="")
		}
		x
	}
	#fun(savtab$scenario[4])
	savtab$scenario<-sapply(savtab$scenario,fun)
	savtab$scenario<-paste("AAA",savtab$scenario,"BBB")
	names(savtab)[1]<-"Scenario"

	knittab<-knitr::kable(
		savtab,
		format="latex",
  		caption=paste("Descriptive statistics for the data simulated under scenarios with", blocknam,". Statistics are mean values for the number of patients under
		tretament and control (trt, ctr), the number of events in either group (ev trt, ev ctr), the number of control group patients who switched (switch), the
		number of patients for whom the time do death was censored by a random censoring event (rcens), the number of simulated data sets out of 10,000 in which the aimed for number
		of events was not achieved within the maximum study duration (ev<aim), the average maximum follow up time in years (FU) and the true hazard ratio conditional on the covariates
		X and W at baseline (HR). For each scenario, 10,000 data sets were simulated.",sep=""),
  		label=paste(names(tab_list)[i],"_descriptive",sep="")
	)
	knittab<-gsub("AAA","\\\\makecell{",knittab)
	knittab<-gsub("MMM","\\\\\\\\",knittab)
	knittab<-gsub("BBB","}",knittab)
	knittab<-gsub("<","$<$",knittab)
	knittab<-gsub(">","$>$",knittab)


	#https://tex.stackexchange.com/questions/318872/automatic-line-breaks-in-a-table

	cat(paste("%",blocknam,"\n"),file=filename,append=i>1)
	cat(knittab,file=filename,append=TRUE)
	cat("\n\n",file=filename,append=TRUE)
}



###
#######common layout:

# global color scale
options(
  ggplot2.discrete.colour = function(...) scale_color_brewer(type="Qualitative", palette = "Set1", ...),
  ggplot2.discrete.fill   = function(...) scale_fill_brewer(type="Qualitative",palette = "Set1", ...)
)

# set global theme
theme_set(
  theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
)

#axis.text=element_text(size=12)
# wrapper to set common graphics device settings
save_plot <- function(gg, filename){
  ggsave(
    filename=filename,
    plot = gg,
    scale=1,
    width=7,
    height=6,
    units="in",
    dpi=600
  )
}

##########
#0.013
#to long

#tab<-tab_list[[1]]$Bias
#tab

#SE<-tab_small_n_H0$SE_bias

make_long<-function(tab) {
	set<-(which(names(tab)=="statistic")+1):(which(names(tab)=="block_nam")-1)
	tab_m<-tab[,set]
	nam<-names(tab_m)
	data.frame(y=unlist(tab_m),Method=rep(nam,each=dim(tab_m)[1]),Scenario=tab$scenario_name)
}

#names(tab_list[[1]])
#"Bias",           "SE_bias"        "Rejection_rate" "CI_coverage"    "CI_width"       "RMSE"

pos <- position_jitterdodge(dodge.width = 0.5, jitter.width = 0)

#setwd("Z:\\EMA_causal_16May_node2\\CI.RCT.Sim\\results4\\Grafiken")


#library(writexl)
#save(tab_list,file="tab_list_21May.RData")

#save(tab_list,file="tab_list_21May.RData")
#save(list=ls(),file="make_tab_list_21May.RData")



farben<-c(
	cens="red",
	ipw="#ff33dd",#"pink",
	ipw_IPCW="pink",
      rpsftm  ="blue",
	rpsftm_rc  ="skyblue",
	rpsftm_rc_IPCW="#87CEED",
      tse = "green2"   ,
	tse_rc ="lightgreen"   ,
	tse_rc_IPCW = "#90EE92",
	gformula = "orange4"     ,
	gformula_IPCW ="orange"
)

#rgb_vals <- col2rgb("lightgreen")
#hex <- rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], maxColorValue = 255)
#hex

#plot(1:length(farben),rep(0,length(farben)),pch=19,col=farben)



#modified from color brewer set 1


#Bias plots
i<-1
farbskala<- scale_colour_manual(values = farben)

titles<-c(
small_n_H1="Small sample size, H1",
small_n_H0="Small sample size, H0",
large_n_H1="Large sample size, H1",
large_n_H0="Large sample size, H0",
extra_cens_H1="Small sample size, H1",
extra_cens_H0="Small sample size, H0"
)



#names(tab_list[[i]])
#ident<-function(x) x
#plot_info<-data.frame(


stat_names=c("Bias","Rejection_rate","CI_coverage","CI_width","RMSE")
ylab_names=c("Relative Bias of HR","Rejection rate","CI coverage","CI width (logHR scale)","RMSE of log HR")

legend_names=c(
 cens="Censor",
ipw ="IPW",
ipw_IPCW="IPW+IPCW",
 rpsftm ="RPSFTM",
rpsftm_rc ="RPSFTM-recensor",
rpsftm_IPCW ="RPSFTM+IPCW",
tse="TSE",
 tse_rc="TSE-recensor",
 tse_IPCW ="TSE+IPCW",
gformula="g-formula",
 gformula_IPCW="g-formula+IPCW")


YLIM<-ylim(0.66,1.5)
REFLINE<- geom_hline(yintercept = 1)
YLAB<-"Relative Bias"

#for(j in 1:5) source("/home/robin/EMA_causal_16May_node2/CI.RCT.Sim/scripts/plot_script_to_be_sourced.R")
for(j in 1:5) source("Z:\\EMA_causal_16May_node2\\CI.RCT.Sim\\scripts\\plot_script_to_be_sourced_1-1.R")


stat<-stat_names[2]
j<-1
j<-2
j<-3
j<-4
j<-5

###############################
#############################

for(j in 1:length(stat_names)) {
  stat<-stat_names[j]
  plot_list<-list()
  levels<-list()



  for(i in 1:length(tab_list)) {
	#plot_data<-make_long(tab_list[[i]]$Bias)
	plot_data<-make_long(tab_list[[i]][[stat]])

	plot_data$Method<-factor(plot_data$Method,levels=names(farben))
	levels[[i]]<-levels(droplevels(plot_data$Method))


	TITEL<-titles[names(tab_list)[i]]
	YLAB<-ylab_names[j]
	if(stat=="Bias") plot_data$y<-exp(plot_data$y)

	plot_rel_bias <- plot_data |>
	    ggplot(aes(x = Scenario, colour = Method, group = Method)) + ylab(YLAB)+
	    #ggplot(aes(x = Scenario, colour = Method, group = Method, shape=Method)) + ylab("Relative Bias")+
	    theme(axis.text.x = element_text(angle =45, hjust = 1,size=7)) +
	    farbskala +
		#scale_colour_brewer(palette = "Set1") +
	    #scale_colour_manual(values = farben) +
	    aes(y = y)  +
	     #geom_line() +
	    geom_point(position = pos) +
	    #geom_errorbar(aes(ymin = y - SE, ymax = y + SE), width = 0.5, position = pos) +
	   #geom_hline(yintercept = REFLINE) + ylim(ylow,yup)+
		scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "),width = 40))+
  theme(
    axis.text.x = element_text(lineheight = 0.6)
  )
	if(stat=="Bias") plot_rel_bias<-plot_rel_bias+geom_hline(yintercept = 1)
	#if(stat=="Rejection_rate" & grepl("H0",TITEL)) plot_rel_bias<-plot_rel_bias+geom_hline(yintercept = 0.025)
	if(stat=="Rejection_rate" & grepl("H0",TITEL)) plot_rel_bias<-plot_rel_bias+geom_hline(yintercept = 0.025)+geom_hline(yintercept = 0.025+sqrt(0.025*(1-0.025)/10000)*qnorm(0.975),linetype = 2)


	plot_list[[i]]<-plot_rel_bias+theme(legend.position = "none")+
		theme(plot.title = element_text(size = 12))+ggtitle(TITEL)
	plot_list
	#save_plot(plot_rel_bias, paste(names(tab_list)[[i]],"_rel_bias.pdf",sep=""))
}


  #pdf("trt_switch_bias.pdf",width=12,height=12)
	filename<-paste("trt_switch_",stat,".pdf",sep="")
  pdf(filename,width=12,height=12)
	plot.new()
	wrap_plots(plot_list[1:4], ncol = 2) + theme(plot.margin = margin(5, 5, 60, 5))
	par(mar=c(0,5,0,5))
	plot.window(new=TRUE,xlim=c(0,1),ylim=c(0,1))
	LEV<-levels[[1]]
	legend("bottom",title="Method",legend=LEV,col=farben[LEV],pch=19,horiz=TRUE,bty="n")
  dev.off()

}
#################

#install.packages("staplr")
#library(staplr)
#remove_pages(
#  rmpages = 1,
#  input_filepath = filename,
#  output_filepath = filename
#)









(p1 | p2 | p3) /
(p4 | p5 | p6) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

library(patchwork)
a<-list(1,2,3)
a[1:2]
plot_list[[2]]
#plot_list[[1]]<-plot_list[[1]]+ theme(legend.position = "bottom")
wrap_plots(plot_list[1:4]) +
  plot_layout(ncol = 2)
wrap_plots(plot_list[1:4]) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

wrap_plots(plots) +
  plot_layout(ncol = 2)

plot.window(new=TRUE,xlim=c(0,1),ylim=c(0,1))
legend("center",legend="A",lwd=1,col=1)
wrap_plots(plots, ncol = 3) +
  plot_layout(
    heights = c(1, 0.8)   # top row, bottom row smaller or larger space control
  )







+
  plot_layout(
    guides = "collect"
  ) &
  theme(legend.position = "bottom")











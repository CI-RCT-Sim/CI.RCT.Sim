save_4by4<-FALSE
stat<-stat_names[j]
plot_list_leg<-list()
  plot_list<-list()
  levels<-list()



  for(i in 1:length(tab_list)) {

	#plot_data<-make_long(tab_list[[i]]$Bias)
	plot_data<-make_long(tab_list[[i]][[stat]])

	SE_tab<-make_long(tab_list[[i]][["SE_bias"]])
	bias_tab<-make_long(tab_list[[i]][["Bias"]])
	SElogHR<-SE_tab$y
	VarTrueVal<-4/(40*10000)
	SEbias<-sqrt(SElogHR^2 + VarTrueVal)
	#SE_Q<-exp(bias_tab$y)*SE_tab$y*qnorm(0.975)
	plot_data$SE_Q<-exp(bias_tab$y)*SEbias*qnorm(0.975)



	#low<-exp(bias_tab$y-SE_tab$y*qnorm(0.975))
	#up<-exp(bias_tab$y+SE_tab$y*qnorm(0.975))


	plot_data$Method<-factor(plot_data$Method,levels=names(farben))
	plot_data$Method<-droplevels(plot_data$Method)
	levels[[i]]<-levels(plot_data$Method)

	plot_data$Scenario<-factor(plot_data$Scenario,levels=unique(plot_data$Scenario))
	#unique(plot_data$Method)
	#table(plot_data$Method)

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
	if(stat=="Bias") plot_rel_bias<-plot_rel_bias+geom_hline(yintercept = 1)+geom_errorbar(aes(ymin = y - SE_Q, ymax = y + SE_Q), width = 0.5, position = pos)
	#if(stat=="Bias") plot_rel_bias<-plot_rel_bias+geom_hline(yintercept = 1)+geom_errorbar(aes(ymin = low, ymax = up), width = 0.5, position = pos)

	#if(stat=="Rejection_rate" & grepl("H0",TITEL)) plot_rel_bias<-plot_rel_bias+geom_hline(yintercept = 0.025)
	if(stat=="Rejection_rate" & grepl("H0",TITEL)) plot_rel_bias<-plot_rel_bias+geom_hline(yintercept = 0.025)+geom_hline(yintercept = 0.025+sqrt(0.025*(1-0.025)/10000)*qnorm(0.975),linetype = 2)

	if(stat=="CI_coverage") {
		plot_rel_bias<-plot_rel_bias+geom_hline(yintercept = 0.95)+geom_hline(yintercept = 0.95+c(-1,1)*sqrt(0.95*(1-0.95)/10000)*qnorm(0.975),linetype = 2)
	}

	plot_list_leg[[i]]<-plot_rel_bias+
		theme(plot.title = element_text(size = 10))+ggtitle(TITEL)
	plot_list[[i]]<-plot_rel_bias+theme(legend.position = "none")+
		theme(plot.title = element_text(size = 12))+ggtitle(TITEL)
	#plot_list
	#save_plot(plot_rel_bias, paste(names(tab_list)[[i]],"_rel_bias.pdf",sep=""))
}


  #pdf("trt_switch_bias.pdf",width=12,height=12)

filename_single<-paste("trt_switch_single_",stat,".pdf",sep="")
#plot_list_leg[[6]]
save_plot(plot_list_leg,file=filename_single)

if(save_4by4) {

	filename<-paste("trt_switch_",stat,".pdf",sep="")
  pdf(filename,width=12,height=12)
	plot.new()
	wrap_plots(plot_list[1:4], ncol = 2) + theme(plot.margin = margin(5, 5, 60, 5))
	par(mar=c(0,5,0,5))
	plot.window(new=TRUE,xlim=c(0,1),ylim=c(0,1))
	LEV<-levels[[1]]
	legend("bottom",title="Method",legend=legend_names[LEV],col=farben[LEV],pch=19,bty="n",ncol=4)
  dev.off()




remove_pages(
  rmpages = 1,
  input_filepath = filename,
  output_filepath = filename
)


	filename_IPCW<-paste("trt_switch_IPCW_",stat,".pdf",sep="")


 pdf(filename_IPCW,width=12,height=8)
	plot.new()
	wrap_plots(plot_list[5:6], ncol = 2) + theme(plot.margin = margin(5, 5, 60, 5))
	par(mar=c(0,5,0,5))
	plot.window(new=TRUE,xlim=c(0,1),ylim=c(0,1))
	LEV<-levels[[5]]
	legnam<-legend_names[LEV]
	legnam<-c(legnam[1],"",legnam[2:9])
	frb<-farben[LEV]
	frb<-c(frb[1],"white",frb[2:9])
	legend("bottom",title="Method",legend=legnam,col=frb,pch=19,bty="n",ncol=5)
  dev.off()

remove_pages(
  rmpages = 1,
  input_filepath = filename_IPCW,
  output_filepath = filename_IPCW
)


}


###############################################################################
# Source code for Section B.4
# Results on Tuning Parameter Selection
###############################################################################



source("./Source/parser.R")


set.seed(1001)

model.no = 3
n = 200
SNR = 2


lambda = seq(0,90,1) * pi / 180



### Example model generation

model.name = paste0("./Sim_data/model",model.no,".R")
source(model.name)
model = model.generate( n = n, U = NULL, SNR = SNR)



### (Load the saved model)


# save.name1 = paste0("./Sim_result/Sim_generic_data_(SecB.4)/m",k,".SNR", SNR ,".data.Rdata")
# load(file = save.name1)



### Estimate signal ranks


r.list = model$r.list

### Estimation
result = estimation.fun(model = model,  lambda = lambda, r.list = r.list)




### Distances from the true structure

dis.test = rep(0,length(lambda))
for(i in 1:length(lambda)){
  dis.test[i] = distance.struct(result$structure.list[[i]]$structure,  model$true.struct)
}





r.test.df = data.frame(lambda = lambda*180/pi, 
                       risk = result$r.test.list, 
                       distance = dis.test,
                       distance.hat = result$diff.list)


### Plot


library(ggplot2)


plot.name = paste0("  Case ",model.no,"  n = ",n,  "  SNR = ", SNR)

scale.df = (max(r.test.df$risk) - min(r.test.df$risk)) /
  (max(c(r.test.df$distance, r.test.df$distance.hat)) - min(c(r.test.df$distance, r.test.df$distance.hat)))

ggplot(r.test.df, aes(x=lambda, y=risk, linetype = "risk")) +
  geom_line() +
  geom_line(data = r.test.df, aes(x = lambda, y = (distance - min(c(distance,distance.hat ))) * scale.df + min(risk),
                                  linetype = "dist") ) +
  geom_line(data = r.test.df, aes(x = lambda, y = (distance.hat - min(c(distance,distance.hat ))) * scale.df + min(risk),
                                  linetype = "dist.h") ) +
  scale_y_continuous(
    
    sec.axis = sec_axis(trans = ~scales::rescale(., to = range(c(r.test.df$distance,r.test.df$distance.hat ))),
                        breaks = function(values) {scales::pretty_breaks(n=5)(values)},
                        name = "dist")) +
  scale_linetype_manual(values = c("risk" = "solid", "dist" = "dashed", "dist.h" = "dotted")) +
  scale_color_manual(values = c("risk" = "black", "dist" = "grey", "dist.h" = "black")) +
  theme_bw()+
  labs(title = plot.name) +
  theme(legend.position = "none")



### Save the results


# save.name1 = paste0("./Sim_result/Sim_generic_data_(SecB.4)/m",model.no,".SNR", SNR ,".data.Rdata")
# save.name2 = paste0("./Sim_result/Sim_generic_data_(SecB.4)/m",model.no,".SNR", SNR ,".result.Rdata")

# save(model, file = save.name1)
# save(result, file = save.name2)



######################################################################
# Draw Figure B.6
# => See "./Sim_codes/Sim_generic_data_plot_(SecB.4).R"
######################################################################


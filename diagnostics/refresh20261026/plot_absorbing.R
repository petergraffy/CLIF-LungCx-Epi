suppressPackageStartupMessages(library(ggplot2))
root <- 'output/run_revision_gaps_20261026'
x <- read.csv(file.path(root, 'sensitivities/absorbing_v2/coordinator/full/prototypes.csv'))
x <- subset(x, domain == 'support')
states <- c('ROOM_AIR','LOW_O2','HFNC','CPAP','NIPPV','IMV','OTHER','UNKNOWN','LIVE_DISCHARGE','HOSPICE','DEATH')
colors <- c('#B3CDE3','#377EB8','#4DAF4A','#984EA3','#6A3D9A','#E6550D','#D8B365','#999999','#009E73','#CC79A7','#222222')
names(colors) <- states
labels <- c('Room air','Low-flow oxygen','HFNC','CPAP','NIPPV','IMV','Other','Unknown','Live hospital discharge','Hospice discharge','Death')
x$state <- factor(x$state, levels=states)
groups <- sort(unique(x$phenotype))
x$group <- factor(x$phenotype, levels=groups, labels=sub('^P','D',groups))
stopifnot(all(is.finite(x$probability)), all(x$probability >= 0),
          max(abs(aggregate(probability~phenotype+bin,x,sum)$probability-1)) < 1e-8)
p <- ggplot(x,aes(bin*6+3,probability,fill=state)) +
  geom_col(width=6,position=position_stack(reverse=TRUE)) +
  facet_wrap(~group,ncol=2) +
  scale_fill_manual(values=colors,breaks=states,labels=labels,drop=FALSE) +
  scale_x_continuous(breaks=c(0,24,48,72),limits=c(0,72),expand=c(0,0)) +
  scale_y_continuous(labels=scales::percent,expand=c(0,0)) +
  labs(title='Respiratory trajectories including early terminal events',
       subtitle='Separate five-group absorbing-state sensitivity | UCMC, NU and Rush',
       x='Hours after first ICU arrival',y='Reconstructed prototype state probability',fill=NULL,
       caption=paste('Death, hospice and live hospital discharge are absorbing states in the underlying hourly sequences.',
                     'Six-hour modal bins and interval reconstruction smooth the aggregate templates; small bands do not establish observed events.',
                     'D groups are descriptive and are not the primary P groups. No outcome regression uses these outcome-containing exposures.',sep='\n')) +
  theme_minimal(base_size=12) +
  theme(legend.position='bottom',plot.caption=element_text(hjust=0,size=9),panel.spacing=grid::unit(1.2,'lines')) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE))
out <- file.path(root,'reports/figures/absorbing_terminal_state_profiles')
ggsave(paste0(out,'.png'),p,width=12,height=10,dpi=180)
ggsave(paste0(out,'.pdf'),p,width=12,height=10)

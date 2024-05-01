library(flowCore)
library(flowAI)
library(openCyto)
library(ggcyto)
library(tidyverse)
library(RColorBrewer)

# Load data as flowSet
fs <- read.flowSet(path = "sample_data", pattern = ".fcs", alter.names = T)
pData(fs)

# View available channels
colnames(fs)

# Raw data processing:
# 1) Automatic QC using flowAI
# 2) Logical transformation for each sample
fs_clean <- flow_auto_qc(fs)
trans <- estimateLogicle(fs_clean[[1]], # Use the unstained sample here
                         colnames(fs_clean[, 3:5])) # Only transform fluorophore channels
fs_clean_trans <- transform(fs_clean, 
                            trans)

# Size Gating to remove debris
# Try a number of different gates to fit your data
g.debris <- polygonGate(filterId = "no_debris",
                        "FSC.A" = c(6e4, 25e4, 25e4, 6e4),
                        "SSC.A" = c(1.5e4, 1.5e4, 25e4, 25e4)) # define gate

# Plot and visualize gating on unstained sample
ggcyto(fs_clean_trans[[1]], aes(x = FSC.A, y = SSC.A), subset = "root") +
  geom_hex(bins = 128) + 
  geom_gate(g.debris) +
  theme_minimal()

# Create a gating set from above boundaries and apply to all samples
gs <- GatingSet(fs_clean_trans) # create a GatingSet
gs_pop_add(gs, g.debris, parent = "root")
recompute(gs)

# Visualize size gating applied to ALL samples
ggcyto(gs,
       aes(x = FSC.A, y = SSC.A),
       subset = "root") + 
  geom_hex(bins = 256) + 
  geom_gate(g.debris, size = 1) +
  ggcyto_par_set(limits = "instrument") +
  facet_wrap(~name, ncol = 3) + 
  scale_fill_gradientn(colours = brewer.pal(n = 8, name = "PiYG"), trans = "sqrt") +
  theme_classic() +
  theme(text = element_text(size = 10))

# Fluorophore gating
# Try an initial gate to begin - can be modified
g.apc <- rectangleGate(filterId = "APC positive",
                       "APC.A" = c(2.4, Inf)) # set gate

# Visualize gate applied to all samples
ggcyto(gs, 
       aes(x = APC.A), 
       subset = "no_debris") +
  geom_density(fill = "magenta") +
  geom_gate(g.apc) +
  geom_stats() +
  ggcyto_par_set(limits = "instrument") +
  facet_wrap(~name, ncol = 2) + 
  theme_classic()

# Add fluorophore gate
gs_pop_add(gs, g.apc, parent = "no_debris") 
recompute(gs) 

# Track gating stats
gs_pop_get_count_with_meta(gs)

# 2D scatter to visualize multiple channels
ggcyto(gs,
       aes(x = FITC.A, y = APC.A),
       subset = "no_debris") + 
  geom_hex(bins = 128) + 
  geom_gate(quadGate(`FITC.A` = 2.5, `APC.A` = 2.4), colour = "black") +
  facet_wrap(~name, ncol = 2) +
  geom_stats() + 
  theme_classic() +
  scale_fill_viridis_c(option = "turbo")

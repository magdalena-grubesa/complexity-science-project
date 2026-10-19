# Multi-symptom alleviating interventions in a depression-anxiety network
   # based on Li et al. (2024), following the logic of NIRA (Lunansky et al., 2022).
   # Main analyses: total dose of 2 SD split over the targeted symptoms
   #   (single symptoms, all pairs, stepwise random, stepwise most influential).
   # Additional analyses: 2 SD per symptom, and a total dose of 4 SD.

############# 1. NETWORK (based on Li et al., 2024) ####################
library(qgraph)
library(ggplot2)

                  

# based on Li et al. (2024)
dep <- c("PHQ1", "PHQ2", "PHQ3", "PHQ4", "PHQ5", "PHQ6", "PHQ7", "PHQ8", "PHQ9")
anx <- c("GAD1", "GAD2", "GAD3", "GAD4", "GAD5", "GAD6", "GAD7")
nodes <- c(dep, anx)
n_nodes <- length(nodes)
is_dep <- nodes %in% dep

labels_long <- c(
  PHQ1 = "Anhedonia", PHQ2 = "Sad mood", PHQ3 = "Sleep", PHQ4 = "Energy",
  PHQ5 = "Appetite", PHQ6 = "Guilt", PHQ7 = "Concentration",
  PHQ8 = "Motor", PHQ9 = "Suicide",
  GAD1 = "Nervousness", GAD2 = "Uncontrollable worry", GAD3 = "Excessive worry",
  GAD4 = "Trouble relaxing", GAD5 = "Restlessness", GAD6 = "Irritability",
  GAD7 = "Afraid")

### EDGE WEIGHTS###
W <- matrix(0.1, n_nodes, n_nodes, dimnames = list(nodes, nodes))   # weak links everywhere
diag(W) <- 0 #symptom cannot have connectivity with itself
# add_edge(a, b, w): sets the connection between symptoms a and b to strength w.
# A connection goes both ways, so it is written in both cells of the table
# (row a / column b AND row b / column a). "<<-" changes the real W, not a temporary copy.
add_edge <- function(a, b, w) { W[a, b] <<- w; W[b, a] <<- w } 

# strong links (thickest lines in Li et al.(2024))
add_edge("PHQ2", "PHQ9", 1.2); add_edge("PHQ2", "PHQ1", 1.2)
add_edge("GAD2", "GAD3", 1.2); add_edge("GAD1", "GAD2", 1.2)
add_edge("PHQ7", "PHQ8", 1.2)

# medium links within depression / within anxiety
add_edge("PHQ1", "PHQ4", 0.7); add_edge("PHQ1", "PHQ7", 0.7)
add_edge("PHQ3", "PHQ4", 0.7); add_edge("PHQ2", "PHQ6", 0.7)
add_edge("PHQ6", "PHQ9", 0.7); add_edge("GAD1", "GAD3", 0.7)
add_edge("GAD1", "GAD4", 0.7); add_edge("GAD2", "GAD4", 0.7)

# bridge links between anxiety and depression (Li et al. (2024), section 3.2)
add_edge("GAD6", "PHQ4", 1.0); add_edge("GAD5", "PHQ8", 1.0)
add_edge("PHQ6", "GAD3", 1.0)
add_edge("GAD4", "PHQ6", 0.5); add_edge("GAD5", "PHQ6", 0.5)
add_edge("GAD6", "PHQ6", 0.5); add_edge("GAD4", "PHQ4", 0.5)

### THRESHOLDS###
# Not published by Li et al. (2024), so symptoms are grouped by how common they
# were in their sample (Table 1, item mean M): more common = easier to activate.
# More negative = harder to switch on. Exact values are our choice;
# the grouping comes from Li et al.'s Table 1 (M).

tau <- setNames(rep(NA, n_nodes), nodes)
tau[c("PHQ6", "GAD6", "GAD3", "PHQ4", "GAD1")] <- -1.6  # very common  (M >= 0.60): guilt, irritability, excessive worry, energy, nervousness
tau[c("PHQ1", "PHQ2", "GAD4", "GAD7")]         <- -2.0  # common       (M 0.50-0.59): anhedonia, sad mood, trouble relaxing, afraid
tau[c("PHQ3", "PHQ5", "PHQ7", "GAD2", "PHQ8")] <- -2.4  # less common  (M 0.40-0.49): sleep, appetite, concentration, uncontrollable worry, motor
tau[c("GAD5", "PHQ9")]                         <- -2.8  # rare         (M < 0.40): restlessness, suicide

tau          # check: all 16 should have a value (no NA)

# Figure: the network
qgraph(W, layout = "spring",
       groups = list(Depression = which(is_dep), Anxiety = which(!is_dep)),
       color = c("blue", "red"), labels = nodes, nodeNames = labels_long,
       legend.cex = 0.4, minimum = 0.1,      # hide the weak 0.1 links
       title = "Network based on Li et al. (2024)")

total_dose <- 2                      # total intervention size in SD (as in Li et al.)

# Function: alleviating intervention. The total dose is split over the targets:
# 1 symptom -> 2 SD, 2 symptoms -> 1 SD each
intervene <- function(thresholds, targets) {
  thresholds[targets] <- thresholds[targets] - (total_dose / length(targets)) * sd(thresholds)
  thresholds
}

############# 2. FUNCTIONS AND UNTREATED NETWORK ####################
# Function: EXACT mean number of active symptoms (no randomness, no seed needed)
# With 16 symptoms there are only 2^16  possible symptom patterns, so we
# calculate the exact average instead of simulating individuals
all_patterns    <- as.matrix(expand.grid(rep(list(0:1), n_nodes)))   # every possible pattern
n_active        <- rowSums(all_patterns)                             # number of symptoms in each pattern
connection_part <- rowSums((all_patterns %*% W) * all_patterns) / 2  # effect of the connections per pattern
exact_mean <- function(thresholds) {
  logp <- as.vector(all_patterns %*% thresholds) + connection_part   # how likely each pattern is (log scale)
  p <- exp(logp - max(logp))                                         # turn into probabilities
  sum(p * n_active) / sum(p)                                         # average number of active symptoms
}

  
original_mean   <- exact_mean(tau)
original_mean

# function: compare one intervention with the original 
compare_with_original <- function(targets) {
  mean_after <- exact_mean(intervene(tau, targets))
  data.frame(
    targets      = paste(targets, collapse = " + "),
    names        = paste(labels_long[targets], collapse = " + "),
    n_symptoms   = length(targets),
    mean_after   = mean_after,
    difference   = mean_after - original_mean)     # negative = improvement
}

############# 3. SINGLE SYMPTOMS (2 SD each, as in Li et al.) ####################
single_results <- data.frame()
for (s in nodes) single_results <- rbind(single_results, compare_with_original(s))
single_results <- single_results[order(single_results$difference), ]
rownames(single_results) <- NULL

single_results


############# 4. ALL 120 PAIRS (total 2 SD = 1 SD each) ####################

all_pairs <- combn(nodes, 2)         # 2 x 120 matrix: every pair of symptoms
ncol(all_pairs)                      # = 120

pair_results <- data.frame()
for (i in 1:ncol(all_pairs)) {
  pair_results <- rbind(pair_results, compare_with_original(all_pairs[, i]))
}
pair_results <- pair_results[order(pair_results$difference), ]   # best to worst pair
rownames(pair_results) <- NULL

head(pair_results, 10)               # 10 best pairs: most negative difference 
tail(pair_results, 10)               # 10 worst pairs

# Figure: the 10 best and 10 worst pairs, with names
top_bottom <- rbind(head(pair_results, 10), tail(pair_results, 10))
top_bottom$group <- rep(c("10 best", "10 worst"), each = 10)
ggplot(top_bottom, aes(x = reorder(names, -difference), y = difference, fill = group)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("10 best" = "lightgreen", "10 worst" = "purple")) +
  labs(x = NULL, y = "Change in mean sum score", fill = NULL,
       title = "Best and worst symptom pairs") +
  theme_minimal()

############# 5. SUMMARY: SINGLE SYMPTOMS VS PAIRS ####################
summary_table <- data.frame(
  intervention    = c("1 symptom (2 SD)", "2 symptoms (1 SD each)"),
  number_tested   = c(nrow(single_results), nrow(pair_results)),
  mean_difference = c(mean(single_results$difference), mean(pair_results$difference)),
  best            = c(min(single_results$difference), min(pair_results$difference)),
  worst           = c(max(single_results$difference), max(pair_results$difference)))

summary_table

# Figure: effect of every pair, sorted from best to worst
pair_results$rank <- 1:nrow(pair_results)               # 1 = best pair
ggplot(pair_results, aes(x = rank, y = difference)) +
  geom_col(fill = "blue") +
  geom_hline(yintercept = min(single_results$difference),
             linetype = "dashed", colour = "red") +     # best single symptom
  labs(x = "Pair (sorted from best to worst)",
       y = "Change in mean sum score",
       title = "All 120 two-symptom interventions (1 SD each)",
       subtitle = "Red dashed line = best single symptom (2 SD)") +
  theme_minimal()

# How many pairs to perform better than the best single symptom?
best_single <- min(single_results$difference)
sum(pair_results$difference < best_single)


############# 6. STEPWISE: RANDOM SYMPTOMS (total 2 SD) ####################
# Same total dose (2 SD) split over k symptoms: k = 1 -> 2 SD, k = 4 -> 0.5 SD each, ...
# For each k we test up to 'n_sets' different random symptom sets
# (for k = 1, 15 and 16 there are only 16, 16 and 1 possible sets, so all are used).
# Takes a few minutes.
set.seed(2)
n_sets <- 100                                   # random symptom sets 

stepwise_results <- data.frame()

for (k in 1:n_nodes) {
  n_possible <- choose(n_nodes, k)             # how many different sets of k symptoms exist
  if (n_possible <= 16) {
    sets <- combn(nodes, k, simplify = FALSE)  # few possible sets: use them all
  } else {
    sets <- replicate(n_sets, sample(nodes, k), simplify = FALSE)   # many: pick random ones
  }
  for (targets in sets) {
    stepwise_results <- rbind(stepwise_results, compare_with_original(targets))
  }}


# Summary per number of symptoms
stepwise_summary <- do.call(rbind, lapply(split(stepwise_results, stepwise_results$n_symptoms),
                                          function(d) data.frame(n_symptoms      = d$n_symptoms[1],
                                                                 sets_tested     = nrow(d),
                                                                 mean_difference = mean(d$difference),
                                                                 best            = min(d$difference),
                                                                 worst           = max(d$difference),
                                                                 spread_sd       = sd(d$difference))))
rownames(stepwise_summary) <- NULL

# Figure: stepwise, random symptoms
ggplot(stepwise_summary, aes(x = n_symptoms, y = mean_difference)) +
  geom_ribbon(aes(ymin = best, ymax = worst), fill = "grey", alpha = 0.2) +  # band = best to worst
  geom_line(colour = "red", linewidth = 1) +
  geom_point(colour = "red", size = 2) +
  scale_x_continuous(breaks = 1:n_nodes) +
  labs(x = "Number of symptoms treated (same total dose of 2 SD, split equally)",
       y = "Change in mean sum score",
       title = "Stepwise: spreading the same treatment over more symptoms",
       subtitle = "Line = average over random symptom sets; band = best to worst set") +
  theme_minimal()



############# 7. STEPWISE: MOST INFLUENTIAL SYMPTOMS (total 2 SD) ####################
# Symptoms are added in order of influence (from section 3: single_results, best first):
# k = 1 -> most influential symptom, k = 2 -> top 2,..k = 16 -> all symptoms
# Same total dose (2 SD) split over the k symptoms, as in the random stepwise

ranking <- single_results$targets                  # symptom codes, most influential first
ranking

influential_results <- data.frame()
for (k in 1:n_nodes) {
  top_k <- ranking[1:k]                            # the k most influential symptoms
  influential_results <- rbind(influential_results, compare_with_original(top_k))
}
influential_results                                # which symptoms were added, and the effect

# Combine with the random stepwise results for one figure
compare_steps <- rbind(
  data.frame(n_symptoms = stepwise_summary$n_symptoms,
             difference = stepwise_summary$mean_difference,
             selection  = "random symptoms (average)"),
  data.frame(n_symptoms = influential_results$n_symptoms,
             difference = influential_results$difference,
             selection  = "most influential symptoms"))

# Figure: most influential vs random symptoms
ggplot(compare_steps, aes(x = n_symptoms, y = difference, colour = selection)) +
  geom_ribbon(data = stepwise_summary, inherit.aes = FALSE,
              aes(x = n_symptoms, ymin = best, ymax = worst),
              fill = "grey", alpha = 0.2) +        # band = best to worst random set
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_colour_manual(values = c("most influential symptoms" = "darkblue",
                                 "random symptoms (average)" = "red")) +
  scale_x_continuous(breaks = 1:n_nodes) +
  labs(x = "Number of symptoms treated (same total dose of 2 SD, split equally)",
       y = "Change in mean sum score",
       colour = NULL,
       title = "Stepwise: most influential vs random symptoms",
       subtitle = "Grey band = best to worst random set") +
  theme_minimal() + theme(legend.position = "bottom")

############# 8. ADDITIONAL: PAIRS WITH 2 SD PER SYMPTOM (4 SD in total) ####################
# Same analysis as before, but every targeted symptom gets the full dose of treatment (2 SD)
# So a pair receives twice as much treatment as a single symptom

intervene_split <- intervene                       
intervene <- function(thresholds, targets) {       # full dose: 2 SD for EACH targeted symptom
  thresholds[targets] <- thresholds[targets] - total_dose * sd(thresholds)
  thresholds
}

# all 120 pairs, each symptom 2 SD
pair_results_full <- data.frame()
for (i in 1:ncol(all_pairs)) {
  pair_results_full <- rbind(pair_results_full, compare_with_original(all_pairs[, i]))
}
pair_results_full <- pair_results_full[order(pair_results_full$difference), ]
rownames(pair_results_full) <- NULL

# summary: single symptom (2 SD) vs pair (2 SD each)
summary_full <- data.frame(
  intervention    = c("1 symptom (2 SD)", "2 symptoms (2 SD each)"),
  mean_difference = c(mean(single_results$difference), mean(pair_results_full$difference)),
  best            = c(min(single_results$difference),  min(pair_results_full$difference)),
  worst           = c(max(single_results$difference),  max(pair_results_full$difference)))
summary_full

# how many pairs (with double the dose) are still WORSE than the best single symptom?
sum(pair_results_full$difference > min(single_results$difference)) 

# Figure: all pairs (2 SD each), with best single symptom as red line
pair_results_full$rank <- 1:nrow(pair_results_full)
ggplot(pair_results_full, aes(x = rank, y = difference)) +
  geom_col(fill = "blue") +
  geom_hline(yintercept = min(single_results$difference),
             linetype = "dashed", colour = "red") +
  labs(x = "Pair (sorted from best to worst)",
       y = "Change in mean sum score",
       title = "All 120 pairs (2 SD each = 4 SD total)",
       subtitle = "Red dashed line = best single symptom (2 SD total)") +
  theme_minimal()

intervene <- intervene_split                       # switch back to the split version

############# 9. ADDITIONAL: STEPWISE WITH 2 SD PER SYMPTOM ####################
# Each targeted symptom gets the full 2 SD, so treating k symptoms = k x 2 SD in total
# Our question: how many symptoms must a treatment hit, and where does the improvement level off?

intervene_split <- intervene                       # store current (split) version
intervene <- function(thresholds, targets) {       # full dose: 2 SD for EACH targeted symptom
  thresholds[targets] <- thresholds[targets] - total_dose * sd(thresholds)
  thresholds
}
n_sets <- 100                                       # random symptom sets per k

set.seed(2)
stepwise_results_full <- data.frame()
for (k in 1:n_nodes) {
  n_possible <- choose(n_nodes, k)                 # how many different sets of k symptoms exist
  if (n_possible <= 16) {
    sets <- combn(nodes, k, simplify = FALSE)      # few possible sets: use them all
  } else {
    sets <- replicate(n_sets, sample(nodes, k), simplify = FALSE)   # many: pick random ones
  }
  for (targets in sets) {
    stepwise_results_full <- rbind(stepwise_results_full, compare_with_original(targets))
  }                       
}

# Summary per number of symptoms
stepwise_summary_full <- do.call(rbind, lapply(split(stepwise_results_full, stepwise_results_full$n_symptoms),
                                               function(d) data.frame(n_symptoms      = d$n_symptoms[1],
                                                                      sets_tested     = nrow(d),
                                                                      mean_difference = mean(d$difference),
                                                                      best            = min(d$difference),
                                                                      worst           = max(d$difference))))
rownames(stepwise_summary_full) <- NULL
stepwise_summary_full$percent_change <- 100 * stepwise_summary_full$mean_difference / original_mean
stepwise_summary_full$extra_gain     <- c(NA, -diff(stepwise_summary_full$mean_difference))  # what each extra symptom adds
stepwise_summary_full

# Figure: stepwise, 2 SD per symptom
ggplot(stepwise_summary_full, aes(x = n_symptoms, y = mean_difference)) +
  geom_ribbon(aes(ymin = best, ymax = worst), fill = "grey", alpha = 0.3) +   # band = best to worst set
  geom_line(colour = "red", linewidth = 1) +
  geom_point(colour = "red", size = 2) +       
  scale_x_continuous(breaks = 1:n_nodes) +
  labs(x = "Number of symptoms treated (each symptom -2 SD)",
       y = "Change in mean sum score",
       title = "Stepwise: broader treatment (2 SD per symptom)",
       subtitle = "Line = average over random sets; band = best to worst") +
       theme_minimal()

intervene <- intervene_split                       # switch back to the split version again


intervene  <- intervene_split        # back to the split version
total_dose <- 4                      # total dose of 4 SD from here on


############# 10. ADDITIONAL: STEPWISE WITH A TOTAL DOSE OF 4 SD ####################
set.seed(2)
stepwise_results_4sd <- data.frame()
for (k in 1:n_nodes) {
  n_possible <- choose(n_nodes, k)
  if (n_possible <= 16) {
    sets <- combn(nodes, k, simplify = FALSE)
  } else {
    sets <- replicate(n_sets, sample(nodes, k), simplify = FALSE)
  }
  for (targets in sets) {
    stepwise_results_4sd <- rbind(stepwise_results_4sd, compare_with_original(targets))
  }
}

stepwise_summary_4sd <- do.call(rbind, lapply(split(stepwise_results_4sd, stepwise_results_4sd$n_symptoms),
  function(d) data.frame(n_symptoms      = d$n_symptoms[1],
                         sets_tested     = nrow(d),
                         mean_difference = mean(d$difference),
                         best            = min(d$difference),
                         worst           = max(d$difference))))
rownames(stepwise_summary_4sd) <- NULL
stepwise_summary_4sd

ggplot(stepwise_summary_4sd, aes(x = n_symptoms, y = mean_difference)) +
  geom_ribbon(aes(ymin = best, ymax = worst), fill = "grey", alpha = 0.2) +
  geom_line(colour = "red", linewidth = 1) +
  geom_point(colour = "red", size = 2) +
  scale_x_continuous(breaks = 1:n_nodes) +
  labs(x = "Number of symptoms treated (same total dose of 4 SD, split equally)",
       y = "Change in mean sum score",
       title = "Stepwise: spreading 4 SD over more symptoms",
       subtitle = "Line = average over random symptom sets; band = best to worst set") +
  theme_minimal()

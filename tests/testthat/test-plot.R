################################################################################
# plot_ordination unit tests
################################################################################
library("phyloseq"); library("testthat"); library("ggplot2")
data("GlobalPatterns")
# Subset to small dataset for quicker testing
GP <- prune_taxa(tail(names(sort(taxa_sums(GlobalPatterns))), 50), GlobalPatterns)

# Pretend GP doesn't have sample_data or tax_table
GP.tax <- tax_table(GP)
GP.sd  <- sample_data(GP)
GP.tr  <- phy_tree(GP)
# GP <- phyloseq(otu_table(GP), GP.tr)
GP.otu <- otu_table(GP)

# Try ordination
GP.ord <- ordinate(GP.otu, "DCA")

# test_that encapsulation makes it difficult to fully test the formula / get()
# step in the formula conversion, but that deprecated workaround is
# still included and will hopefully bridge the gap for users switching
# from previous formula-first use-cases, where the left-hand side of
# the formula specified the phyloseq-data
#test_that("plot_ordination: formula-first should give a deprecation warning", {
# expect_warning(GP.ord.cap <- ordinate(GP~SampleType, "CAP"))
# expect_warning(GP.ord.cca <- ordinate(GP~SampleType, "CCA"))
# expect_warning(GP.ord.rda <- ordinate(GP~SampleType, "RDA"))
# # But it still works.
# expect_is(GP.ord.cap, "capscale")
# expect_is(GP.ord.cca, "cca")
# expect_is(GP.ord.rda, "rda")


test_that("plot_ordination: Naked otu_table results in warning, but no error", {
  expect_is(GP.ord, "decorana")
	# samples-only
	expect_that(plot_ordination(GP.otu, GP.ord, "samples"), gives_warning())
	# species. 
	expect_that(plot_ordination(GP.otu, GP.ord, "species"), gives_warning())
	# split
	expect_that(plot_ordination(GP.otu, GP.ord, "split"), gives_warning())
	# biplot
	expect_that(plot_ordination(GP.otu, GP.ord, "biplot"), gives_warning())
})

# Create (merged) phyloseq-class GP, and run comparisons
test_that("all 4 plot_ordination type options result in valid ggplot2 object", {
	GP <- merge_phyloseq(GP.otu, GP.tr)
	# Print. Don't want the render directive to have an error,
  # even while the ggplot object is created.
	expect_is(print(plot_ordination(GP, GP.ord, "samples")), "gg")
	expect_is(print(plot_ordination(GP, GP.ord, "species")), "gg")
	expect_is(print(plot_ordination(GP, GP.ord, "split")), "gg")
	expect_is(print(plot_ordination(GP, GP.ord, "biplot")), "gg")
	# Don't print. Test that result is ggplot-class
	expect_is(plot_ordination(GP, GP.ord, "samples"), "ggplot")
	expect_is(plot_ordination(GP, GP.ord, "species"), "ggplot")
	expect_is(plot_ordination(GP, GP.ord, "split"), "ggplot")
	expect_is(plot_ordination(GP, GP.ord, "biplot"), "ggplot")
})

test_that("plot_ordination: The justDF=TRUE option returns a data.frame", {
  # Make GP a phyloseq object with only tree (no ordination co-variables to plot)
	GP <- merge_phyloseq(GP.otu, GP.tr)
	expect_that(df0 <- plot_ordination(GP, GP.ord, "species", justDF=TRUE), is_a("data.frame"))	
	expect_that(df1 <- plot_ordination(GP, GP.ord, "samples", justDF=TRUE), is_a("data.frame"))	
	expect_that(df2 <- plot_ordination(GP, GP.ord, "split", justDF=TRUE), is_a("data.frame"))		
	expect_that(df3 <- plot_ordination(GP, GP.ord, "biplot", justDF=TRUE), is_a("data.frame"))	
	# split and biplot data.frames should be same.	
	expect_that(df2, is_identical_to(df3))
})

test_that("plot_ordination: When variables are present or not, color SampleType", {
	p1 <- plot_ordination(GP, GP.ord, "samples", color="SampleType")	
	expect_that(p2<-plot_ordination(GP, GP.ord, "species", color="SampleType"), gives_warning())
	p3 <- plot_ordination(GP, GP.ord, "split", color="SampleType")	
	p4 <- plot_ordination(GP, GP.ord, "biplot", color="SampleType")
	# ggplot-class tests
	expect_is(p1, "ggplot")
	expect_is(p2, "ggplot")
	expect_is(p3, "ggplot")
	expect_is(p4, "ggplot")
	expect_is(print(p1), "gg")
	expect_is(print(p2), "gg")
	expect_is(print(p3), "gg")
	expect_is(print(p4), "gg")
})


test_that("plot_ordination: When variables are present or not, shape SamplyType", {
	# GP <- merge_phyloseq(GP.otu, GP.tr, GP.sd, GP.tax)
	# Pair down samples to just five sampleTypes, for shape plotting.
	GP <- subset_samples(GP, SampleType %in% c("Feces", "Freshwater", "Ocean",
                                             "Tongue", "Sediment (estuary)"))
	# Some legend issues here that need tidying...
	p1 <- plot_ordination(GP, GP.ord, "samples", shape="SampleType")	
	expect_warning(p2 <- plot_ordination(GP, GP.ord, "species", shape="SampleType"))
	p3 <- plot_ordination(GP, GP.ord, "split", shape="SampleType")	
	p4 <- plot_ordination(GP, GP.ord, "biplot", shape="SampleType")	
	# ggplot-class tests
	expect_is(p1, "ggplot")
	expect_is(p2, "ggplot")
	expect_is(p3, "ggplot")
	expect_is(p4, "ggplot")
	expect_is(print(p1), "gg")
	expect_is(print(p2), "gg")
	expect_is(print(p3), "gg")
	expect_is(print(p4), "gg")	
})

test_that("plot_ordination: When variables are present or not, label SamplyType", {
	p1 <- plot_ordination(GP, GP.ord, "samples", label="SampleType")	
	expect_warning(p2 <- plot_ordination(GP, GP.ord, "species", label="SampleType"))
	p3 <- plot_ordination(GP, GP.ord, "split", label="SampleType")	
	p4 <- plot_ordination(GP, GP.ord, "biplot", label="SampleType")
	# ggplot-class tests
	expect_is(p1, "ggplot")
	expect_is(p2, "ggplot")
	expect_is(p3, "ggplot")
	expect_is(p4, "ggplot")
	expect_is(print(p1), "gg")
	expect_is(print(p2), "gg")
	expect_is(print(p3), "gg")
	expect_is(print(p4), "gg")	
})

test_that("plot_ordination: Continuous variables still mapped, uses added dummy variable", {
	# Add the fake continuous variable
	sample_data(GP)$OMEGA3_FA_CONC <- sample(1:100, nsamples(GP)) 
	expect_is(p1 <- plot_ordination(GP, GP.ord, "samples", color="OMEGA3_FA_CONC"), "ggplot")
  # Continuous variable cannot be mapped to shape. This is a ggplot object,
  # but will throw error when 'printed'
	p2 <- plot_ordination(GP, GP.ord, "samples", shape="OMEGA3_FA_CONC")
	expect_is(p2, "ggplot")
  expect_error(print(p2))
  # A `label` can be mapped to continuous var. It is coerced to char and printed.
	expect_is(p3 <- plot_ordination(GP, GP.ord, "samples", label="OMEGA3_FA_CONC"), "ggplot")
	expect_that(print(p1), is_a("gg"))
	#expect_that(print(p2), is_a("gg"))
	expect_that(print(p3), is_a("gg"))			
})

test_that("plot_ordination: Some additional formats and warnings.", {
  GP.ord.cca = ordinate(GP, "CCA")
  GP.ord.mdsbray = ordinate(GP, "MDS", "bray")
  expect_is(p1 <- plot_ordination(GP, GP.ord.mdsbray, type="TaXa", color="Phylum", title="p1"), "ggplot")
  expect_is(p2 <- plot_ordination(GP, GP.ord.mdsbray, type="SpLit", color="Phylum", title="p2"), "ggplot")
  expect_warning(p3 <- plot_ordination(GP, GP.ord.mdsbray, type="SamplE", color="Kingdom", title="p3"))
  expect_is(p3, "ggplot")
  expect_is(p4 <- plot_ordination(GP, GP.ord.cca, type="TaXa", color="Kingdom", title="p4"), "ggplot")
  expect_is(p5 <- plot_ordination(GP, GP.ord.cca, type="samPle", color="SampleType", title="p5"), "ggplot")
  expect_is(p6 <- plot_ordination(GP, GP.ord.cca, type="biplot", color="SampleType", title="p6"), "ggplot")
  expect_is(p7 <- plot_ordination(GP, GP.ord.cca, type="biplot",
                         label="X.SampleID", color="SampleType", title="p7"), "ggplot")
  expect_is(p7b <- plot_ordination(GP, GP.ord.cca, type="biplot",
                          label="X.SampleID", color=NULL, title="p7b"), "ggplot")
  expect_is(p7c <- plot_ordination(GP, GP.ord.cca, type="biplot",
                          label="Phylum", color=NULL, title="p7c"), "ggplot")
  expect_is(p7d <- plot_ordination(GP, GP.ord.cca, type="biplot",
                          label="Phylum", color="SampleType", title="p7d"), "ggplot")
  expect_is(p8 <- plot_ordination(GP, GP.ord.cca, type="scree",
                         label="X.SampleID", color="SampleType", title="p8"), "ggplot")
  expect_is(p9 <- plot_ordination(GP, GP.ord.cca, type=" sPlit __ ",
                         label="Phylum", color="SampleType", title="p8"), "ggplot")
  expect_that(print(p1), is_a("gg"))
  expect_that(print(p2), is_a("gg"))
  expect_that(print(p3), is_a("gg"))
  expect_that(print(p4), is_a("gg"))
  expect_that(print(p5), is_a("gg"))
  expect_that(print(p6), is_a("gg"))
  expect_that(print(p7), is_a("gg"))
  expect_that(print(p7b), is_a("gg"))
  expect_that(print(p7c), is_a("gg"))
  expect_that(print(p7d), is_a("gg"))
  expect_that(print(p8), is_a("gg"))
  expect_that(print(p9), is_a("gg"))
  # A few more related to new `wascores` support as default backup coordinates
  xnames = tapply(taxa_sums(GlobalPatterns), tax_table(GlobalPatterns)[, "Phylum"], sum)
  xnames <- names(sort(xnames, decreasing = TRUE))[1:5]
  GP = prune_taxa(taxa_sums(GlobalPatterns) > 1E4, GlobalPatterns)
  GP <- prune_taxa(tax_table(GP)[, "Phylum"] %in% xnames, GP)
  x = ordinate(GP, method = "MDS", distance = "unifrac", weighted=TRUE)
  y = ordinate(GP, method = "CCA")
  z = ordinate(GP, method = "CAP", "unifrac", ~SampleType)
  z1 = ordinate(GP, method = "CAP", "bray", ~SampleType)
  # Try a bunch more with splits and biplots
  expect_is(p11 <- plot_ordination(GP, x, type = "biplot", color="Phylum"), "ggplot")
  expect_is(print(p11), "gg")
  expect_is(p12 <- plot_ordination(GP, x, type = "biplot", color="SampleType", shape="Phylum"), "ggplot")
  expect_is(print(p12), "gg")
  expect_is(p13 <- plot_ordination(GP, x, type = "split", color="SampleType", shape="Phylum"), "ggplot")
  expect_is(print(p13), "gg")
  expect_is(p14 <- plot_ordination(GP, x, type = "split", color="Phylum"), "ggplot")
  expect_is(print(p14), "gg")
  expect_is(p15 <- plot_ordination(GP, y, type = "biplot", color="Phylum"), "ggplot")
  expect_is(print(p15), "gg")
  expect_is(p16 <- plot_ordination(GP, y, type = "species", color="Phylum"), "ggplot")
  expect_is(print(p16), "gg")
  expect_is(p17 <- plot_ordination(GP, z, type = "biplot", color="Phylum"), "ggplot")
  expect_is(print(p17), "gg")
  expect_is(p18 <- plot_ordination(GP, z, type = "biplot", color="SampleType", shape="Phylum"), "ggplot")
  expect_is(print(p18), "gg")
  expect_is(p19 <- plot_ordination(GP, z1, type = "biplot", color="Phylum"), "ggplot")
  expect_is(print(p19), "gg")
})

test_that("plot_ordination: CAP method", {
  # Works with a named formula argument
  GP.ord.cap1 = ordinate(GP, method="CAP", distance="bray", formula=~SampleType)
  expect_is(GP.ord.cap1, "capscale")
  # Works without naming the formula argument
  GP.ord.cap2 = ordinate(GP, method="CAP", distance="bray", ~SampleType)
  expect_is(GP.ord.cap2, "capscale")
  expect_equivalent(GP.ord.cap1, GP.ord.cap2)  
  # Works with precomputed distance matrix
  Dist = distance(GP, "bray", type="samples")
  GP.ord.cap3 = ordinate(physeq=GP, method="CAP", distance=Dist, ~SampleType)
  expect_is(GP.ord.cap3, "capscale")
  # Can't expect equivalent b/c pre-computed distance
  # won't carryover any species/taxa scores.
  #expect_equivalent(GP.ord.cap1, GP.ord.cap3)
  GP.ord.cap = ordinate(GP, method="CAP", distance="bray", formula=~SampleType)
  expect_is(GP.ord.cap, "capscale")
  expect_is(p4 <- plot_ordination(GP, GP.ord.cap, type="TaXa",
                                  color="Phylum", title="p4"), "ggplot")
  expect_is(p5 <- plot_ordination(GP, GP.ord.cap, type="samPle", 
                                  color="SampleType", title="p5"), "ggplot")
  expect_is(p6 <- plot_ordination(GP, GP.ord.cap, type="biplot",
                                  color="SampleType", title="p6"), "ggplot")
  expect_is(p7 <- plot_ordination(GP, GP.ord.cap, type="biplot", label="X.SampleID",
                                  color="SampleType", title="p7"), "ggplot")
  expect_is(p7b <- plot_ordination(GP, GP.ord.cap, type="biplot", label="X.SampleID",
                                   color=NULL, title="p7b"), "ggplot")
  expect_is(p7c <- plot_ordination(GP, GP.ord.cap, type="biplot",
                                   label="Phylum", color=NULL, title="p7c"), "ggplot")
  expect_is(p7d <- plot_ordination(GP, GP.ord.cap, type="biplot", label="Phylum",
                                   color="SampleType", title="p7d"), "ggplot")
  expect_is(p8 <- plot_ordination(GP, GP.ord.cap, type="scree", label="X.SampleID",
                                  color="SampleType", title="p8"), "ggplot")
  expect_is(p9 <- plot_ordination(GP, GP.ord.cap, type=" sPlit __ ", label="Phylum",
                                  color="SampleType", title="p8"), "ggplot")
  expect_that(print(p4), is_a("gg"))
  expect_that(print(p5), is_a("gg"))
  expect_that(print(p6), is_a("gg"))
  expect_that(print(p7), is_a("gg"))
  expect_that(print(p7b), is_a("gg"))
  expect_that(print(p7c), is_a("gg"))
  expect_that(print(p7d), is_a("gg"))
  expect_that(print(p8), is_a("gg"))
  expect_that(print(p9), is_a("gg"))
})

# Constrained CCA / RDA
test_that("plot_ordination: CCA, RDA method", {
  # Constrained RDA and CCA both work.
  GP.ord.cca = ordinate(GP, "CCA", NULL, formula=~SampleType)
  expect_is(GP.ord.cca, "cca")
  GP.ord.rda = ordinate(GP, "RDA", NULL, formula=~SampleType)
  expect_is(GP.ord.rda, "rda")
  # Test plotting CCA
  expect_is(p4 <- plot_ordination(GP, GP.ord.cca, type="TaXa",
                                  color="Phylum", title="p4"), "ggplot")
  expect_is(p5 <- plot_ordination(GP, GP.ord.cca, type="samPle", 
                                  color="SampleType", title="p5"), "ggplot")
  expect_is(p6 <- plot_ordination(GP, GP.ord.cca, type="biplot",
                                  color="SampleType", title="p6"), "ggplot")
  expect_is(p7 <- plot_ordination(GP, GP.ord.cca, type="biplot", label="X.SampleID",
                                  color="SampleType", title="p7"), "ggplot")
  expect_is(p7b <- plot_ordination(GP, GP.ord.cca, type="biplot", label="X.SampleID",
                                   color=NULL, title="p7b"), "ggplot")
  expect_is(p7c <- plot_ordination(GP, GP.ord.cca, type="biplot",
                                   label="Phylum", color=NULL, title="p7c"), "ggplot")
  expect_is(p7d <- plot_ordination(GP, GP.ord.cca, type="biplot", label="Phylum",
                                   color="SampleType", title="p7d"), "ggplot")
  expect_is(p8 <- plot_ordination(GP, GP.ord.cca, type="scree", label="X.SampleID",
                                  color="SampleType", title="p8"), "ggplot")
  expect_is(p9 <- plot_ordination(GP, GP.ord.cca, type=" sPlit __ ", label="Phylum",
                                  color="SampleType", title="p8"), "ggplot")
  expect_that(print(p4), is_a("gg"))
  expect_that(print(p5), is_a("gg"))
  expect_that(print(p6), is_a("gg"))
  expect_that(print(p7), is_a("gg"))
  expect_that(print(p7b), is_a("gg"))
  expect_that(print(p7c), is_a("gg"))
  expect_that(print(p7d), is_a("gg"))
  expect_that(print(p8), is_a("gg"))
  expect_that(print(p9), is_a("gg"))
  # Repeat test-plotting RDA
  expect_is(p4 <- plot_ordination(GP, GP.ord.rda, type="TaXa",
                                  color="Phylum", title="p4"), "ggplot")
  expect_is(p5 <- plot_ordination(GP, GP.ord.rda, type="samPle", 
                                  color="SampleType", title="p5"), "ggplot")
  expect_is(p6 <- plot_ordination(GP, GP.ord.rda, type="biplot",
                                  color="SampleType", title="p6"), "ggplot")
  expect_is(p7 <- plot_ordination(GP, GP.ord.rda, type="biplot", label="X.SampleID",
                                  color="SampleType", title="p7"), "ggplot")
  expect_is(p7b <- plot_ordination(GP, GP.ord.rda, type="biplot", label="X.SampleID",
                                   color=NULL, title="p7b"), "ggplot")
  expect_is(p7c <- plot_ordination(GP, GP.ord.rda, type="biplot",
                                   label="Phylum", color=NULL, title="p7c"), "ggplot")
  expect_is(p7d <- plot_ordination(GP, GP.ord.rda, type="biplot", label="Phylum",
                                   color="SampleType", title="p7d"), "ggplot")
  expect_is(p8 <- plot_ordination(GP, GP.ord.rda, type="scree", label="X.SampleID",
                                  color="SampleType", title="p8"), "ggplot")
  expect_is(p9 <- plot_ordination(GP, GP.ord.rda, type=" sPlit __ ", label="Phylum",
                                  color="SampleType", title="p8"), "ggplot")
  expect_that(print(p4), is_a("gg"))
  expect_that(print(p5), is_a("gg"))
  expect_that(print(p6), is_a("gg"))
  expect_that(print(p7), is_a("gg"))
  expect_that(print(p7b), is_a("gg"))
  expect_that(print(p7c), is_a("gg"))
  expect_that(print(p7d), is_a("gg"))
  expect_that(print(p8), is_a("gg"))
  expect_that(print(p9), is_a("gg"))
})
################################################################################
# Other plot function tests...
################################################################################
# plot_richness tests
################################################################################
test_that("estimate_richness: test values, classes", {
  data("soilrep")
  data("GlobalPatterns")
  # Default is all available measures
  erdf = estimate_richness(soilrep)
  expect_is(erdf, "data.frame")
  expect_equivalent(nrow(erdf), 56)
  # Contains all expected measures 
  expect_true(all(c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher") %in% colnames(erdf)))
  # and certain standard errors:
  expect_true(all(c("se.chao1", "se.ACE") %in% colnames(erdf)))
  # Test some values. 
  expect_equivalent(erdf$Observed, apply(otu_table(soilrep), 2, function(x){sum(x>0)}))
  expect_equivalent(estimate_richness(GlobalPatterns, measures="Observed")[, 1], 
                    apply(otu_table(GlobalPatterns), 2, function(x){sum(x>0)}))
  # Calculate "manually" the values that should be Chao1, compare with result.
  S_0 = apply(otu_table(soilrep), 2, function(x){sum(x>0)})
  a1  = apply(otu_table(soilrep), 2, function(x){sum(x==1)})
  a2  = apply(otu_table(soilrep), 2, function(x){sum(x==2)})
  S_P = S_0 + a1*(a1-1)/(2*(a2+1))
  expect_equivalent(round(S_P, 4), round(estimate_richness(soilrep, measures="Chao1")[, "Chao1"], 4))
  # Expect a data.frame, even with just one column
  expect_is(estimate_richness(soilrep, measures="Observed"), "data.frame")
  # Specify a few:
  x = estimate_richness(soilrep, measures=c("Observed", "InvSimpson", "Shannon", "Chao1"))
  expect_equivalent(round(x[1:5, "Shannon"], 4), 
                    round(c(6.540578, 6.715170, 6.948412, 7.343088, 6.838917), 4))
})

test_that("plot_richness: Standard plots work", {
  data("soilrep")
  p = plot_richness(soilrep)
  expect_is(p, "ggplot")
  expect_equivalent(levels(p$data$variable),
                    c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
  expect_false(all(is.na(p$data$se)))
  expect_true(any(is.na(p$data$se)))
  p = plot_richness(soilrep, measures=c("Observed", "Chao1"))
  expect_is(p, "ggplot")
  expect_equivalent(levels(p$data$variable), c("Observed", "Chao1"))
})

test_that("plot_richness: sortby argument works correctly", {
  data("soilrep")
  # sortby must be among the `measures`.
  # Should throw warning if not, but still produce a plot.
  expect_warning({p1 <- plot_richness(soilrep, sortby="Treatment")})  
  expect_is(p1, "ggplot")
  # sortby is only relevant if `x` argument is discrete.
  # Should throw warning if not, but still produce a plot.
  # First add dummy numeric sample variable
  sample_data(soilrep)$dummy <- runif(nsamples(soilrep))
  expect_warning({p2 <- plot_richness(soilrep, x="dummy", sortby="Chao1")})
  expect_is(p2, "ggplot")
  # Default `x` is "samples", always discrete
  p3 = plot_richness(soilrep, sortby="Chao1")
  expect_equivalent(levels(p3$data$samples)[1:5],
                    c("a_C137", "a_C145", "a_C126", "a_C156", "a_C139"))
  expect_is(p3, "ggplot")
  # Make sure the discrete aggregation sort gets the order correct as well.
  p4 = plot_richness(soilrep, x="Treatment", sortby="Simpson") 
  expect_equivalent(levels(p4$data$Treatment), c("UC", "WC", "WU", "UU"))
  expect_is(p4, "ggplot")
})

test_that("plot_richness/estimate_richness: fisher.alpha", {
  data("GlobalPatterns")
  data("soilrep")
  p = plot_richness(soilrep, measures="Fisher")
  expect_is(p, "ggplot")
  expect_is(p123123 <- plot_richness(GlobalPatterns, measures="Fisher"), "ggplot")
  expect_equivalent(levels(p123123$data$variable), "Fisher")
})
################################################################################
# Test psmelt properly protects against various name collisions
################################################################################
test_that("psmelt properly protects against various name collisions", {
  data("GlobalPatterns")
  gp.ch = subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
  ps1 = NULL
  gp1 = gp.ch
  # type-1a conflict, Abundance
  sample_data(gp1)$Abundance <- paste0("Sa-", 1:nsamples(gp1))
  expect_warning(ps1 <- psmelt(gp1))
  expect_equal(colnames(ps1)[1:3], c("OTU", "Sample", "Abundance"))
  expect_equal(dim(ps1), c(546L, 18L))
  expect_true("sample_Abundance" %in% colnames(ps1))
  # A different type-1a conflict, OTU
  ps1 = NULL
  gp1 = gp.ch
  sample_data(gp1)$OTU <- paste0("Sa-", 1:nsamples(gp1))
  expect_warning(ps1 <- psmelt(gp1))
  expect_equal(colnames(ps1)[1:3], c("OTU", "Sample", "Abundance"))
  expect_equal(dim(ps1), c(546L, 18L))
  expect_true("sample_OTU" %in% colnames(ps1))
  # A different type-1a conflict, Sample
  ps1 = NULL
  gp1 = gp.ch
  sample_data(gp1)$Sample <- paste0("Sa-", 1:nsamples(gp1))
  expect_warning(ps1 <- psmelt(gp1))
  expect_equal(colnames(ps1)[1:3], c("OTU", "Sample", "Abundance"))
  expect_equal(dim(ps1), c(546L, 18L))
  expect_true("sample_Sample" %in% colnames(ps1))  
  # type-1b conflict. rank_names conflict with special variables
  ps1 = NULL
  gp1 = gp.ch
  tax_table(gp1) <- cbind(tax_table(gp1), Sample=paste0("ta", taxa_names(gp1)))
  expect_warning(ps1 <- psmelt(gp1))
  expect_equal(colnames(ps1)[1:3], c("OTU", "Sample", "Abundance"))
  expect_equal(dim(ps1), c(546L, 18L))
  expect_true("taxa_Sample" %in% colnames(ps1))  
  # type-2 conflict. Variable collision between rank_names and sample_data
  ps1 = NULL
  gp1 = gp.ch
  tax_table(gp1) <- cbind(tax_table(gp1), Primer=paste0("ta", taxa_names(gp1)))
  expect_warning(ps1 <- psmelt(gp1))
  expect_equal(colnames(ps1)[1:3], c("OTU", "Sample", "Abundance"))
  expect_equal(dim(ps1), c(546L, 18L))
  expect_true("sample_Primer" %in% colnames(ps1)) 
  # All conflict types at once.
  ps1 = NULL
  gp1 = gp.ch
  sample_data(gp1)$Abundance <- paste0("Sa-", 1:nsamples(gp1))  
  sample_data(gp1)$OTU <- paste0("Sa-", 1:nsamples(gp1))
  sample_data(gp1)$Sample <- paste0("Sa-", 1:nsamples(gp1))
  tax_table(gp1) <- cbind(tax_table(gp1), Sample=paste0("ta", taxa_names(gp1)))
  tax_table(gp1) <- cbind(tax_table(gp1), Primer=paste0("ta", taxa_names(gp1)))
  expect_warning(ps1 <- psmelt(gp1))
  expect_equal(colnames(ps1)[1:3], c("OTU", "Sample", "Abundance"))
  expect_equal(dim(ps1), c(546L, 22L))
  newvars = c("sample_OTU", "sample_Sample", "sample_Abundance",
              "sample_Primer", "taxa_Sample")
  expect_true(all(newvars %in% colnames(ps1)))   
})
################################################################################
test_that("psmelt correctly handles phyloseq data with NULL components, and OTU tables", {
  data("GlobalPatterns")
  GP = prune_taxa(names(sort(taxa_sums(GlobalPatterns), TRUE)[1:50]), GlobalPatterns)
  # The objects with NULL components
  GPS = phyloseq(otu_table(GP), sample_data(GP), phy_tree(GP))
  GPT = phyloseq(otu_table(GP), tax_table(GP), phy_tree(GP))
  GPTr = phyloseq(otu_table(GP), phy_tree(GP))
  GPN = otu_table(GP)
  # Try psmelt directly. Should be no errors or warnings.
  expect_is((testT <- psmelt(GPT)), "data.frame")
  expect_is((testS <- psmelt(GPS)), "data.frame")
  expect_is((testTr <- psmelt(GPTr)), "data.frame")
  expect_is((testN <- psmelt(GPN)), "data.frame")
  # Test values of the results.
  expect_is(testT$Abundance, "numeric")
  expect_is(testT$OTU, "character")
  expect_is(testT$Sample, "character")
  expect_equivalent(colnames(testT), c("OTU", "Sample", "Abundance", "Kingdom", "Phylum",
                                       "Class", "Order", "Family", "Genus", "Species"))
  expect_equivalent(colnames(testS), c("Sample", "OTU", "Abundance", "X.SampleID", "Primer",
                                       "Final_Barcode", "Barcode_truncated_plus_T", 
                                       "Barcode_full_length", "SampleType", "Description"))
  # Try psmelt via plot function that relies on it
  expect_is(pS <- plot_tree(GPS, color="SampleType"), "ggplot")
  expect_is(pT <- plot_tree(GPT, shape="Kingdom"), "ggplot")
  expect_is(pTr <- plot_tree(GPTr), "ggplot")
  expect_is(pN <- plot_bar(GPN), "ggplot")
  expect_is((prPS<-print(pS)), "gg")
  expect_is((prPT<-print(pT)), "gg")
  expect_is((prPTr<-print(pTr)), "gg")
  expect_is((prPN<-print(pN)), "gg")
})
test_that("psmelt doesn't break when the number of taxa is 1", {
  data(GlobalPatterns)
  # tree removal warning when prune to 1 OTU.
  expect_warning(GP1 <- prune_taxa(taxa_names(GlobalPatterns)[1], GlobalPatterns))
  expect_equal(ntaxa(GP1), 1)
  df <- psmelt(GP1)
  expect_is(df, 'data.frame')
  reqnames = c("OTU", "Sample", "Abundance", "SampleType", "Kingdom", "Phylum")
  expect_true(all(reqnames %in% names(df)))
  expect_equivalent(sum(df$Abundance, na.rm = TRUE), taxa_sums(GP1))
})
################################################################################

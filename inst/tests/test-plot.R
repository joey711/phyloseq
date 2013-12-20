################################################################################
# plot_ordination unit tests
################################################################################
library("phyloseq"); library("testthat"); library("ggplot2")
data("GlobalPatterns")
# Subset to small dataset for quicker testing
GP <- prune_taxa(taxa_sums(GlobalPatterns)>10000, GlobalPatterns)

# Pretend GP doesn't have sample_data or tax_table
GP.tax <- tax_table(GP)
GP.sd  <- sample_data(GP)
GP.tr  <- phy_tree(GP)
# GP <- phyloseq(otu_table(GP), GP.tr)
GP.otu <- otu_table(GP)

# Try ordination
GP.ord <- ordinate(GP.otu, "DCA")

test_that("plot_ordination: Naked otu_table results in warning, but no error", {
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

	# print. Don't want the render directive to have an error, even while the ggplot object is created.
	expect_that(print(plot_ordination(GP, GP.ord, "samples")), is_a("list"))
	expect_that(print(plot_ordination(GP, GP.ord, "species")), is_a("list"))
	expect_that(print(plot_ordination(GP, GP.ord, "split")), is_a("list"))
	expect_that(print(plot_ordination(GP, GP.ord, "biplot")), is_a("list"))

	# don't print. Test that result is ggplot-class
	expect_that(plot_ordination(GP, GP.ord, "samples"), is_a("ggplot"))
	expect_that(plot_ordination(GP, GP.ord, "species"), is_a("ggplot"))
	expect_that(plot_ordination(GP, GP.ord, "split"), is_a("ggplot"))
	expect_that(plot_ordination(GP, GP.ord, "biplot"), is_a("ggplot"))
})

test_that("plot_ordination: The justDF=TRUE option returns a data.frame", {
	GP <- merge_phyloseq(GP.otu, GP.tr)
	expect_that(df0 <- plot_ordination(GP, GP.ord, "species", justDF=TRUE), is_a("data.frame"))	
	expect_that(df1 <- plot_ordination(GP, GP.ord, "samples", justDF=TRUE), is_a("data.frame"))	
	expect_that(df2 <- plot_ordination(GP, GP.ord, "split", justDF=TRUE), is_a("data.frame"))		
	expect_that(df3 <- plot_ordination(GP, GP.ord, "biplot", justDF=TRUE), is_a("data.frame"))	
	# split and biplot data.frames should be same.	
	expect_that(df2, is_identical_to(df3))
})

test_that("plot_ordination: When variables are present or not, color SampleType", {
	# The full-featured version
	# GP <- merge_phyloseq(GP.otu, GP.tr, GP.sd, GP.tax)
	
	p1 <- plot_ordination(GP, GP.ord, "samples", color="SampleType")	
	p2 <- plot_ordination(GP, GP.ord, "species", color="SampleType")	
	p3 <- plot_ordination(GP, GP.ord, "split", color="SampleType")	
	p4 <- plot_ordination(GP, GP.ord, "biplot", color="SampleType")	
			
	# ggplot-class tests
	expect_that(p1, is_a("ggplot"))
	expect_that(p2, is_a("ggplot"))
	expect_that(p3, is_a("ggplot"))
	expect_that(p4, is_a("ggplot"))
		
	expect_that(print(p1), is_a("list"))
	expect_that(print(p2), throws_error())
	expect_that(print(p3), is_a("list"))
	expect_that(print(p4), is_a("list"))		
})


test_that("plot_ordination: When variables are present or not, shape SamplyType", {
	# GP <- merge_phyloseq(GP.otu, GP.tr, GP.sd, GP.tax)
	# Pair down samples to just five sampleTypes, for shape plotting.
	GP <- subset_samples(GP, SampleType %in% c("Feces", "Freshwater", "Ocean", "Tongue", "Sediment (estuary)"))
	
	# Some legend issues here that need tidying...
	p1 <- plot_ordination(GP, GP.ord, "samples", shape="SampleType")	
	p2 <- plot_ordination(GP, GP.ord, "species", shape="SampleType")	
	p3 <- plot_ordination(GP, GP.ord, "split", shape="SampleType")	
	p4 <- plot_ordination(GP, GP.ord, "biplot", shape="SampleType")	
			
	# ggplot-class tests
	expect_that(p1, is_a("ggplot"))
	expect_that(p2, is_a("ggplot"))
	expect_that(p3, is_a("ggplot"))
	expect_that(p4, is_a("ggplot"))
		
	expect_that(print(p1), is_a("list"))
	expect_that(print(p2), throws_error())
	expect_that(print(p3), is_a("list"))
	expect_that(print(p4), is_a("list"))	
})

test_that("plot_ordination: When variables are present or not, label SamplyType", {
	# GP <- merge_phyloseq(GP.otu, GP.tr, GP.sd, GP.tax)
			
	# ggplot-class tests
	expect_that(p1 <- plot_ordination(GP, GP.ord, "samples", label="SampleType"), is_a("ggplot"))
	expect_that(p2 <- plot_ordination(GP, GP.ord, "species", label="SampleType"), throws_error())
	expect_that(p3 <- plot_ordination(GP, GP.ord, "split", label="SampleType"), is_a("ggplot"))
	expect_that(p4 <- plot_ordination(GP, GP.ord, "biplot", label="SampleType"), is_a("ggplot"))
		
	expect_that(print(p1), is_a("list"))
	# expect_that(print(p2), throws_error())
	expect_that(print(p3), is_a("list"))
	expect_that(print(p4), is_a("list"))
	
})

test_that("plot_ordination: Continuous variables still mapped, uses added dummy variable", {
	# GP <- merge_phyloseq(GP.otu, GP.tr, GP.sd, GP.tax)
	# Add the fake continuous variable
	sample_data(GP)$OMEGA3_FA_CONC <- sample(1:100, nsamples(GP)) 
	
	expect_that(p1 <- plot_ordination(GP, GP.ord, "samples", color="OMEGA3_FA_CONC"), is_a("ggplot"))
	expect_that(p2 <- plot_ordination(GP, GP.ord, "samples", shape="OMEGA3_FA_CONC"), is_a("ggplot"))
	expect_that(p3 <- plot_ordination(GP, GP.ord, "samples", label="OMEGA3_FA_CONC"), is_a("ggplot"))
	
	expect_that(print(p1), is_a("list"))
	expect_that(print(p2), throws_error())
	expect_that(print(p3), is_a("list"))			
})
################################################################################
# Other plot function tests should follow here...
################################################################################

# plot_richness tests
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

test_that("plot_richness/estimate_richness: fisher.alpha", {
  data("GlobalPatterns")
  data("soilrep")
  p = plot_richness(soilrep, measures="Fisher")
  expect_is(p, "ggplot")
  expect_is(p123123 <- plot_richness(GlobalPatterns, measures="Fisher"), "ggplot")
  expect_equivalent(levels(p123123$data$variable), "Fisher")
})

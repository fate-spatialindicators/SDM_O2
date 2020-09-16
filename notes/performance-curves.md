# Various notes on (thermal) performance curves

"For example, thermal biologists have quantified performance curves for species from different latitudes, often concluding that species near the tropics where temperatures are generally high and relatively stable – tend to live closer to their thermal optimum than do species living in temperate latitudes (Deutsch et al., 2008; Morley et al., 2012; Stillman and Somero, 2000), where short-term and seasonal variation is greater."
- @dowd2015

"These conclusions are generally robust to different choices of the function used to describe performance curves, so long as the curve is asymmetrical and left-skewed, including modified BoltzmannArrhenius (Dell et al., 2011), modified beta (Niehaus et al., 2012) (Fig. 1A), Gaussian×Gompertz (Martin and Huey, 2008) and Briere three-parameter (Estay et al., 2014)."
- @dowd2015

"Although previously popular because of their simple parameterization, symmetrical Gaussian performance curves generate dramatically different results and should be avoided for describing performance (Asbury and Angilletta, 2010; Dell et al., 2011; Gilchrist, 1995)."
- @dowd2015

"For example, we have shown that mean thermal tolerance in intertidal limpets appears to be set by rare events of elevated body temperature that are unlikely to occur within any individual’s lifespan (Denny and Dowd, 2012). Similarly, others have concluded that acute thermal tolerance may be more relevant to survival in natural environments than responses to chronic exposures (Angilletta, 2009)."
- @dowd2015

"However, numerous studies have demonstrated that variation over very small scales can rival or even exceed mean differences observed over much larger scales (e.g. Bartlett and Gates, 1967; Denny et al., 2011; Elvin and Gonor, 1979; Miller et al., 2009; Pincebourde and Woods, 2012; Seabra et al., 2011). For example, in the rocky intertidal zone, the difference in body temperatures between the warmest and coolest mussels over an area of a few square meters (up to 15°C on any given day) rivaled and sometimes greatly exceeded the expected difference in body temperatures along ∼1600 km of the western coastline of North America (Denny et al., 2011)."
- @dowd2015

@angilletta2006 fits 5 curves and discusses comparing with AIC etc.:
- Gaussian
- Quadratic
- Weibull
- Modified Gaussian (Estimates exponent instead of ^2)
- Exponentially modified Gaussian (need to look at more carefully)

@niehaus2012 modified beta function (usually best fit) and Weibull (best in some cases)
"We compared six functions: quadratic, Gaussian, modified Gaussian, exponentially modified Gaussian, Weibull and beta functions (supplementary material Table S1). Three of these functions – the Gaussian, quadratic and Weibull functions – have been used to theoretically or empirically describe thermal reaction norms (Huey and Kingsolver, 1993; Palaima and Spitze, 2004). The remaining functions were chosen because their complex structure should provide a better fit to non-linear data."

which is awfully similar to:

@angilletta2006 "Three of these functions—the Gaussian, Quadratic, and Weibull functions—have been used to theoretically or empirically describe thermal performance curves (e.g., see Huey and Kingsolver, 1993; Huang and Yang, 1995; Palaima and Spitze, 2004). The remaining functions—the modiﬁed Gaussian and the exponentially modiﬁed Gaussian functions—were chosen because their complex structure should provide a better ﬁt to nonlinear data."

@martin2008: Gaussian multiplied by a Gompertz function. They include an index to calculate the level of asymmetry.

@luhring2016: Fits Lactin-2 curves. Look viable for us.

@quinn2017: Fits 33(!) thermal performance curves (for arthropods). See the giant table.

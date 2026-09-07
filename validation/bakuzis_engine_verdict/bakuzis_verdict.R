suppressMessages(library(data.table))
d <- fread('/fs/scratch/PUOM0008/crsfaaron/wt-engine/engine_bakuzis_out/engine_bakuzis_100yr.csv')
# QMD from BA(ft2/ac) and TPA: QMD=sqrt(BA*576/(pi*TPA))
d[, qmd := ifelse(tpa>0, sqrt(atba*576/(pi*tpa)), NA_real_)]
d[, sdi := ifelse(tpa>0 & qmd>0, tpa*(qmd/10)^1.605, 0)]
# Engine SDIMAX (from .out SDI MAX table): site species approx
sdimax_map <- c('Spruce-Fir'=506,'Northern-Hardwood'=560,'Pine'=529,'Oak-Pine'=560)
d[, sdimax := sdimax_map[species_group]]
cat('=== per-group max SDI reached over 100yr (should not exceed engine SDIMAX) ===\n')
tt <- d[, .(sdi_max_reached=round(max(sdi,na.rm=T)), sdimax_engine=sdimax[1]), by=species_group]
tt[, overshoot := round(sdi_max_reached/sdimax_engine,2)]
print(tt)

cat('\n=== SELF-THINNING: does SDI approach but not exceed SDIMAX? (flag if >1.05) ===\n')
st <- d[, .(peak_ratio=max(sdi,na.rm=T)/sdimax[1]), by=.(species_group,site_class,density_class)]
cat('n scenarios with SDI>1.05*SDIMAX:', st[peak_ratio>1.05,.N], 'of', st[,.N],'\n')
cat('median peak SDI/SDIMAX:', round(median(st$peak_ratio),2),' max:',round(max(st$peak_ratio),2),'\n')

cat('\n=== SITE ORDERING: higher site -> taller top height at yr100 & more BA growth ===\n')
so <- d[year==2100, .(topht=mean(attopht), ba=mean(atba)), by=.(species_group,site_class)]
so[, site_ord := factor(site_class, levels=c('Low','Medium','High'))]
so <- so[order(species_group, site_ord)]
print(so)
# monotone check per group
mono_site <- so[, .(topht_mono = all(diff(topht[order(site_ord)])>=0)), by=species_group]
cat('\nTop-height increases with site class (per group):\n'); print(mono_site)

cat('\n=== EICHHORN: top height tracks site x age (independent of density) ===\n')
# For a given site+species, top height at fixed age should be ~invariant to initial density
eich <- d[year %in% c(2050,2100), .(topht_cv = sd(attopht)/mean(attopht)), by=.(species_group,site_class,year)]
cat('mean within-(site,age) topht CV across densities:', round(mean(eich$topht_cv,na.rm=T),3),
    ' (small = Eichhorn holds)\n')

cat('\n=== MONOTONE TOP HEIGHT: does engine top height ever DECLINE within a run? ===\n')
d <- d[order(scenario, year)]
d[, dtop := attopht - shift(attopht), by=scenario]
declines <- d[!is.na(dtop) & dtop < -0.5]
cat('scenarios with any topht decline >0.5ft:', length(unique(declines$scenario)),'of 36\n')
cat('total declining steps:', nrow(declines),'of', nrow(d[!is.na(dtop)]),'\n')
if(nrow(declines)>0){ cat('worst declines:\n'); print(head(declines[order(dtop),.(scenario,species_group,site_class,year,attopht,dtop)],8)) }

cat('\n=== YIELD-DENSITY / productivity by site (mean BA yr100) ===\n')
print(d[year==2100, .(mean_ba=round(mean(atba),1), mean_topht=round(mean(attopht),1)), by=site_class][order(site_class)])

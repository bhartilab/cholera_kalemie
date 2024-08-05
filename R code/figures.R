fig1 <- function(){
  # Shapefiles are freely available on the GADM website (https://gadm.org/)
  africa_df <- rgdal::readOGR("PATH TO THE FOLDER/afr_g2014_2013_0.shp")

  drc_df <- africa_df[africa_df$ADM0_NAME=="Democratic Republic of the Congo",]

  africa_df <- africa_df %>%
                sf::st_as_sf()


  # Shapefiles are freely available on the GRIP website (https://www.globio.info/download-grip-dataset)
  road_df <- rgdal::readOGR("PATH TO THE FOLDER/GRIP_DRC.shp")
  road_df <- road_df[road_df$GP_RSE=="1",] %>%
              sf::st_as_sf()

  # Shapefiles are freely available on the Landscan website (https://landscan.ornl.gov/)
  pop_raster <- raster::raster("PATH TO THE FOLDER/landscan-global-2013.tif") %>%
          raster::crop(., drc_df) %>%
          raster::mask(., drc_df) %>%
          as(., "SpatialPixelsDataFrame") %>%
          as.data.frame() %>%
          setNames(., c("value", "x", "y"))

  drc_df <- drc_df %>%
              ggplot2::fortify(., region="ADM0_NAME")
  
  # The zs object comes from shapefiles freely available on the GADM website (https://gadm.org/)
  # After modifyin its attributes to modify how the names of the health zones were managed:
  # the name of the health zones is now 'ZS' and the names are in full capital letters
  # We then saved the modified shapefile as a SpatialPolygonsDataFrame in a rda file 'zs.rda'
  # Download the shapefiles and modify the attributes accordingly to run the following chunk
  load("PATH TO THE FOLDER/zs.rda")

  zs_kal_df <- sf::st_as_sf(zs) %>%
            dplyr::filter(ZS %in% c("KALEMIE", "NYEMBA"))

  zs_df <- sf::st_as_sf(zs)

  zs_raster <- raster::raster("PATH TO THE FOLDER/landscan-global-2013.tif") %>%
          raster::crop(., zs[zs$ZS %in% c("KALEMIE", "NYEMBA"),]) %>%
          raster::mask(., zs[zs$ZS %in% c("KALEMIE", "NYEMBA"),]) %>%
          as(., "SpatialPixelsDataFrame") %>%
          as.data.frame() %>%
          setNames(., c("value", "x", "y"))


  # The shapefiles for the lakes are freely available on the Geoserver website (http://geoserver01.uit.tufts.edu/wfs?outputformat=SHAPE-ZIP&request=GetFeature&service=wfs&srsName=EPSG%3A4326&typeName=sde%3AGISPORTAL.GISOWNER01.DRCWATERPOLY09&version=2.0.0)
  lakes <- rgdal::readOGR("PATH TO THE FOLDER/Africa_waterbody.shp"))
  lakes <- lakes[lakes$NAME_OF_WA %in% c("Tanganyika", "Victoria", "Albert", "Edward", "Kivu", "Mweru", "Malawi", "Turkana", "Rukwa (north)", "Rukwa (south)", "Kyoga"),] %>%
              sf::st_as_sf()

  box_zs <- sf::bbox(zs_kal_df)

  # Data coming from the Supplementary Information of the article "Modalities and preferred routes of geographic spread of cholera from endemic areas in eastern Democratic Republic of the Congo" from Kayembe et al (DOI: 10.1371/journal.pone.0263160)
  data <- read.delim("PATH TO FOLDER/surv_drc.txt", header=TRUE) %>%
      setNames(., c("ID", "ZS", "AN", "NUMSEM", "TOTALCAS")) %>%
      dplyr::mutate(ZS=toupper(ZS)) %>%
      dplyr::filter(ZS %in% c("KALEMIE", "NYEMBA") & AN %in% 2002:2014) %>%
      dplyr::mutate(
        datesem=ISOweek::ISOweek2date(paste0(AN, "-W", stringr::str_pad(NUMSEM, 2, side="left", pad = "0"), "-3"))) %>%
      dplyr::group_by(datesem) %>%
      dplyr::summarize(cases=sum(TOTALCAS)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        w=lubridate::isoweek(datesem))

  africa_map <- ggplot2::ggplot() +
                  ggplot2::geom_sf(
                    data=africa_df,
                    ggplot2::aes(geometry=geometry),
                    color="grey",
                    fill="grey") +
                  ggplot2::geom_sf(
                    data=drc_df,
                    ggplot2::aes(geometry=geometry),
                    color="black",
                    fill="white",
                    size=0.3) +
                  ggplot2::xlab("") +
                  ggplot2::ylab("") +
                  ggplot2::theme_void() +
                  ggplot2::theme(
                    panel.border=ggplot2::element_rect(colour="black", fill=NA))

  drc_map_alt <- ggplot2::ggplot() +
              ggplot2::geom_tile(
                data=pop_raster,
                ggplot2::aes(x=x, y=y, fill=log10(value+1))) +
              ggplot2::geom_sf(
                data=zs_df,
                ggplot2::aes(geometry=geometry),
                color="white",
                fill=NA,
                size=0.2) +
              ggplot2::geom_sf(
                data=lakes %>%
                  dplyr::filter(NAME_OF_WA=="Tanganyika"),
                ggplot2::aes(geometry=geometry),
                fill="#16537e",
                color=NA) +
              ggplot2::geom_sf(
                data=road_df,
                ggplot2::aes(geometry=geometry),
                color="springgreen",
                size=0.5) +
              ggplot2::geom_rect(
                ggplot2::aes(xmin=box_zs[1], xmax=box_zs[3], ymin=box_zs[2], ymax=box_zs[4]),
                color="red",
                fill=NA,
                size=1) +
              ggplot2::geom_point(
                ggplot2::aes(x=mean(sp::coordinates(kalemie)[,1]), y=mean(sp::coordinates(kalemie)[,2])),
                size=2,
                color="red") +
              ggplot2::scale_fill_gradientn(
                name="",
                colors = viridis::viridis_pal(begin=0.2, option="B")(36)) +
              ggplot2::theme_void() +
              ggplot2::theme(
                legend.position=c(0.12, 0.2),
                panel.background = ggplot2::element_rect(
                  fill="white",
                  color=NA),
                legend.title=ggplot2::element_text(size=10))

  zs_map <- ggplot2::ggplot() +
              ggplot2::geom_sf(
                data=lakes,
                ggplot2::aes(geometry=geometry),
                color=NA,
                fill="#558ec0") +
              ggplot2::geom_sf(
                data=zs_df,
                ggplot2::aes(geometry=geometry),
                color="darkgrey",
                fill="grey85",
                linewidth=1) +
              ggplot2::geom_tile(
                data=zs_raster,
                ggplot2::aes(x=x, y=y, fill=log10(value+1))) +
              ggplot2::geom_sf(
                data=zs_kal_df,
                ggplot2::aes(geometry=geometry),
                color="red4",
                fill=NA,
                size=1) +
              ggplot2::geom_point(
                ggplot2::aes(x=mean(sp::coordinates(kalemie)[,1]), y=mean(sp::coordinates(kalemie)[,2])),
                shape=21,
                size=6,
                stroke=1.5,
                color="black",
                fill="red") +
              ggplot2::scale_fill_gradientn(
                name="",
                colours=colorRampPalette(colors = c("lightgrey", "red"))(36)) +
              ggplot2::coord_cartesian(
                xlim=box_zs[1,],
                ylim=box_zs[2,]) +
              ggplot2::theme_void() +
              ggplot2::theme(
                legend.position="none",
                panel.background = ggplot2::element_rect(
                  fill="white",
                  color=NA),
                legend.title=ggplot2::element_text(size=10))


  drc_map_alt <- drc_map_alt +
             ggplot2::annotation_custom(
                ggplot2::ggplotGrob(africa_map), 
                  x = 12.5,
                  y = 1, 
                  xmax = 17.5,
                  ymax = 5.5)

  time_series_alt_eps <- ggplot2::ggplot() +
                  ggplot2::geom_rect(
                    ggplot2::aes(
                      xmin=9,
                      xmax=22,
                      ymin=-Inf,
                      ymax=Inf),
                    fill="grey",
                    color=NA) +
                  ggplot2::geom_rect(
                    ggplot2::aes(
                      xmin=36,
                      xmax=53,
                      ymin=-Inf,
                      ymax=Inf),
                    fill="grey",
                    color=NA) +
                  ggplot2::geom_boxplot(
                    data=data,
                    ggplot2::aes(x=w, y=cases, group=w),
                    fill="tomato") +
                  ggplot2::ylab("Weekly suspected cholera cases") +
                  ggplot2::xlab("ISO weeks") +
                  ggplot2::theme_bw()

  ggplot2::ggsave(
    "PATH TO YOUR FOLDER/fig1_man_A.eps",
    drc_map_alt,
    height=4,
    width=4)

  ggplot2::ggsave(
    "PATH TO YOUR FOLDER/fig1_man_B.eps",
    zs_map,
    height=4,
    width=3.4)

  ggplot2::ggsave(
    "PATH TO YOUR FOLDER/fig1_man_C.eps",
    time_series_alt_eps,
    height=4,
    width=8)
  # The subfigures were put together using adobe illustrator to get the final figure available in the manuscript 
}


fig2 <- function(){
  # This chunk extracts the estimates and the credible intervals from calculations done in "other_functions_migrN.R"
  # The results of this chunk are already available in minimal data sets in /data (see README.md)

  # load("PATH TO THE FOLDER WHERE YOU WANT TO SAVE THE ESTIMATES/estimate_c1.rda")
  # estimate <- estimate_c1
  # rm("estimate_c1")
  # fit_df <- lapply(1:length(estimate), function(x){
  #   return(
  #     estimate[[x]]$timeseries %>%
  #     dplyr::mutate(
  #       sample=x,
  #       N=estimate[[x]]$Nt$N[seq(0, 118, by=0.01) %in% 1:118]))
  # }) %>%
  #   do.call("rbind", .) %>%
  #     dplyr::group_by(t) %>%
  #     dplyr::summarize(
  #       N=mean(N),
  #       lb=bayestestR::ci(I, method="HDI")$CI_low,
  #       up=bayestestR::ci(I, method="HDI")$CI_high,
  #       I=mean(I))

  # prop <- lapply(1:length(estimate), function(x){
  #   return(
  #     estimate[[x]]$prop %>%
  #     dplyr::filter(t %in% seq(0, 118, by=1)) %>%
  #     dplyr::mutate(
  #       sample=x))
  # }) %>%
  #   do.call("rbind", .)
  # rm("estimate")

  load("PATH ON YOUR MACHINE/data/fit_df.rda")
  load("PATH ON YOUR MACHINE/data/prop.rda")
  load("PATH ON YOUR MACHINE/data/cases.rda")

  y <- cases
  date_iso <- data.frame(
                year=c(rep(2013, length(47:52)), rep(2014, 52), rep(2015, 53), rep(2016, 7)),
                week=c(47:52, 1:52, 1:53, 1:7)) %>%
                dplyr::mutate(
                  datesem=ISOweek::ISOweek2date(paste0(year, "-W", stringr::str_pad(week, 2, side="left", pad = "0"), "-3")))

  date_iso_l <- rbind(
                  data.frame(
                    year=2013,
                    week=46,
                    datesem=lubridate::dmy("13/11/2013")),
                  date_iso) %>%
                dplyr::mutate(
                  t=0:118)

  season_df <- rbind(
                data.frame(
                  xbeg=lubridate::dmy(
                    c("13/11/2013", "1/3/2014", "1/9/2014", "1/3/2015", "1/9/2015")),
                  xend=lubridate::dmy(
                    c("31/12/2013", "1/5/2014", "31/12/2014", "1/5/2015", "31/12/2015")),
                  compartment="Susceptibles"),
                data.frame(
                  xbeg=lubridate::dmy(
                    c("13/11/2013", "1/3/2014", "1/9/2014", "1/3/2015", "1/9/2015")),
                  xend=lubridate::dmy(
                    c("31/12/2013", "1/5/2014", "31/12/2014", "1/5/2015", "31/12/2015")),
                  compartment="Infected"),
                data.frame(
                  xbeg=lubridate::dmy(
                    c("13/11/2013", "1/3/2014", "1/9/2014", "1/3/2015", "1/9/2015")),
                  xend=lubridate::dmy(
                    c("31/12/2013", "1/5/2014", "31/12/2014", "1/5/2015", "31/12/2015")),
                  compartment="Recovered")
              )

  precision_eps <- ggplot2::ggplot()  +
      ggplot2::geom_rect(
        data=season_df,
        ggplot2::aes(
          xmin=xbeg,
          xmax=xend,
          ymin=-Inf,
          ymax=Inf),
        fill="grey",
        # alpha=0.2,
        color=NA) +
      ggplot2::geom_ribbon(
        data=prop %>%
          dplyr::select(
            t, Susceptibles=s, Infected=i) %>%
          dplyr::mutate(Recovered=1-Susceptibles-Infected) %>%
          tidyr::gather(., compartment, value, Susceptibles:Recovered) %>%
          dplyr::group_by(t, compartment) %>%
          dplyr::summarize(
            lb=bayestestR::ci(100*value, method="HDI")$CI_low,
            up=bayestestR::ci(100*value, method="HDI")$CI_high,
            value=100*mean(value)) %>%
          dplyr::ungroup() %>%
          dplyr::left_join(
            .,
            date_iso_l %>%
              dplyr::select(t, datesem),
            by="t"),
        ggplot2::aes(x=datesem, ymin=lb, ymax=up),
        fill="gray35",
        # alpha=0.2,
        color=NA) +
      ggplot2::geom_line(
        data=prop %>%
          dplyr::select(
            t, Susceptibles=s, Infected=i) %>%
          dplyr::mutate(Recovered=1-Susceptibles-Infected) %>%
          tidyr::gather(., compartment, value, Susceptibles:Recovered) %>%
          dplyr::group_by(t, compartment) %>%
          dplyr::summarize(
            lb=bayestestR::ci(100*value, method="HDI")$CI_low,
            up=bayestestR::ci(100*value, method="HDI")$CI_high,
            value=100*mean(value)) %>%
          dplyr::ungroup() %>%
          dplyr::left_join(
            .,
            date_iso_l %>%
              dplyr::select(t, datesem),
            by="t"),
        ggplot2::aes(x=datesem, y=value),
        color="gray25",
        size=1) +
      viridis::scale_color_viridis(
        name="",
        option="D",
        discrete=TRUE) +
      viridis::scale_fill_viridis(
        name="",
        option="D",
        discrete=TRUE) +
      ggplot2::facet_wrap(
        .~compartment,
        ncol=1,
        scales="free_y") +
      ggplot2::geom_vline(
        xintercept=date_iso$datesem[1],
        color="darkviolet",
        linetype="dashed") +
      ggplot2::geom_vline(
        xintercept=date_iso$datesem[32],
        color="darkviolet",
        linetype="dashed") +
      ggplot2::geom_vline(
        xintercept=date_iso$datesem[35],
        color="darkviolet",
        linetype="dashed") +
      ggplot2::ylab("% of the population") +
      ggplot2::xlab("") +
      ggplot2::scale_x_date(
        labels=scales::label_date_short("%b-%Y"),
        breaks="years",
        limits=range(date_iso_l$datesem)) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.title.x=ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(size = 26),
        axis.title.y=ggplot2::element_text(size = 26),
        axis.text.y = ggplot2::element_text(size = 26),
        strip.text.x = ggplot2::element_text(
          size = 26,
          margin=ggplot2::margin(t=5, b=5)),
        legend.position="none")

  fit_eps <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin=lubridate::dmy(
          c("13/11/2013", "1/3/2014", "1/9/2014", "1/3/2015", "1/9/2015")),
        xmax=lubridate::dmy(
          c("31/12/2013", "1/5/2014", "31/12/2014", "1/5/2015", "31/12/2015")),
        ymin=-Inf,
        ymax=Inf),
      fill="grey",
      color=NA) +
    ggplot2::geom_ribbon(
      data=fit_df %>%
        dplyr::mutate(datesem=date_iso$datesem),
      ggplot2::aes(x=datesem, ymin=lb, ymax=up),
      color=NA,
      fill="gray35") +
    ggplot2::geom_line(
      data=fit_df %>%
        dplyr::mutate(datesem=date_iso$datesem),
      ggplot2::aes(x=datesem, y=I),
      color="gray25",
      size=1) +
    ggplot2::geom_point(
      data=data.frame(
        t=1:118,
        y=y,
        datesem=date_iso$datesem),
      ggplot2::aes(x=datesem, y=y),
      size=1,
      fill=NA,
      shape=21) +
    ggplot2::geom_vline(
      xintercept=date_iso$datesem[1],
      color="darkviolet",
      linetype="dashed") +
    ggplot2::geom_vline(
      xintercept=date_iso$datesem[32],
      color="darkviolet",
      linetype="dashed") +
    ggplot2::geom_vline(
      xintercept=date_iso$datesem[35],
      color="darkviolet",
      linetype="dashed") +
    ggplot2::facet_wrap(
      .~"Model fit",
      ncol=1) +
    ggplot2::xlab("") +
    ggplot2::ylab("Weekly reported\nsuspected cholera cases") +
    ggplot2::scale_x_date(
      labels=scales::label_date_short("%b-%Y"),
      breaks="years",
      limits=range(date_iso_l$datesem)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
       axis.title.x=ggplot2::element_blank(),
       axis.text.x = ggplot2::element_blank(),
       axis.ticks.x = ggplot2::element_blank(),
       strip.text.x = ggplot2::element_text(
        size = 20,
        margin=ggplot2::margin(t=5, b=5)),
       axis.title.y=ggplot2::element_text(size = 26),
       axis.text.y = ggplot2::element_text(size = 26))

  ggplot2::ggsave(
    "PATH TO YOUR FOLDER/fig2.eps",
    cowplot::plot_grid(
      fit_eps,
      precision_eps,
      ncol=1,
      align="v",
      axis="lr",
      rel_heights=c(1,3),
      labels=c("", "")),
    width=10,
    height=8)
  # The final figure available in the manuscript was then improved in adobe illustrator
}

fig3 <- function(){
  # This chunk extracts the estimates and the credible intervals from calculations done in "other_functions.R"
  # The results of this chunk are already available in minimal data sets in /data (see README.md)

  # load("PATH TO THE FOLDER WHERE YOU WANT TO SAVE THE ESTIMATES/estimate_c1.rda")
  # estimate <- estimate_c1
  # rm("estimate_c1")
  # prop <- lapply(1:length(estimate), function(x){
  #   return(
  #     estimate[[x]]$prop %>%
  #     dplyr::filter(t %in% c(0:118)) %>%
  #     dplyr::mutate(
  #       sample=x))
  # }) %>%
  #   do.call("rbind", .)

  load("PATH ON YOUR MACHINE/data/prop.rda")

  int_data <- rbind(
        prop %>%
          dplyr::select(
            t, Susceptibles=s, Infected=i) %>%
          dplyr::mutate(Recovered=1-Susceptibles-Infected) %>%
          tidyr::gather(., compartment, value, Susceptibles:Recovered) %>%
          dplyr::group_by(t, compartment) %>%
          dplyr::summarize(
            value=mean(value)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(type="Intervention as\nit happened"),
        prop %>%
          dplyr::select(
            t, Susceptibles=s_v, Infected=i_v) %>%
          dplyr::mutate(Recovered=1-Susceptibles-Infected) %>%
          tidyr::gather(., compartment, value, Susceptibles:Recovered) %>%
          dplyr::group_by(t, compartment) %>%
          dplyr::summarize(
            value=mean(value)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(type="No vaccination"),
        prop %>%
          dplyr::select(
            t, Susceptibles=s_w, Infected=i_w) %>%
          dplyr::mutate(Recovered=1-Susceptibles-Infected) %>%
          tidyr::gather(., compartment, value, Susceptibles:Recovered) %>%
          dplyr::group_by(t, compartment) %>%
          dplyr::summarize(
            value=mean(value)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(type="No WASH"),
        prop %>%
          dplyr::select(
            t, Susceptibles=s_w_v, Infected=i_w_v) %>%
          dplyr::mutate(Recovered=1-Susceptibles-Infected) %>%
          tidyr::gather(., compartment, value, Susceptibles:Recovered) %>%
          dplyr::group_by(t, compartment) %>%
          dplyr::summarize(
            value=mean(value)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(type="No vaccination\nand no WASH")
        ) %>%
        dplyr::mutate(
          compartment=factor(
            compartment,
            levels=c("Infected", "Recovered", "Susceptibles")),
          type=factor(
            type,
            levels=c("Intervention as\nit happened",
              "No vaccination",
              "No WASH",
              "No vaccination\nand no WASH")))

  # This chunk extracts the estimates and the credible intervals from calculations done in "other_functions.R"
  # The results of this chunk are already available in minimal data sets in /data (see README.md)

  # cum_df <- lapply(1:length(estimate), function(x){
  #           return(
  #             estimate[[x]]$timeseries %>%
  #             dplyr::mutate(
  #               sample=x))
  #         }) %>%
  #         do.call("rbind", .) %>%
  #         dplyr::group_by(sample) %>%
  #         dplyr::summarize(
  #           i=sum(i),
  #           i_v=sum(i_v),
  #           i_w=sum(i_w),
  #           i_w_v=sum(i_w_v)) %>%
  #         dplyr::mutate(
  #           red_v=i_v-i,
  #           red_w=i_w-i,
  #           red_w_v=i_w_v-i) %>%
  #         dplyr::ungroup() %>%
  #         tidyr::gather(., type, value, red_v:red_w_v) %>%
  #         dplyr::mutate(
  #           type=factor(
  #             type,
  #             levels=c("red_v", "red_w", "red_w_v"),
  #             labels=c("No vaccination",
  #               "No WASH",
  #               "No vaccination\nand no WASH")))

  load("PATH ON YOUR MACHINE/data/cum_df.rda")

  date_iso <- data.frame(
                year=c(rep(2013, length(47:52)), rep(2014, 52), rep(2015, 53), rep(2016, 7)),
                week=c(47:52, 1:52, 1:53, 1:7)) %>%
                dplyr::mutate(
                  datesem=ISOweek::ISOweek2date(paste0(year, "-W", stringr::str_pad(week, 2, side="left", pad = "0"), "-3")))

  date_iso_l <- rbind(
                  data.frame(
                    year=2013,
                    week=46,
                    datesem=lubridate::dmy("13/11/2013")),
                  date_iso) %>%
                dplyr::mutate(
                  t=0:118)

  season_df <- rbind(
                data.frame(
                  xbeg=lubridate::dmy(
                    c("13/11/2013", "1/3/2014", "1/9/2014", "1/3/2015", "1/9/2015")),
                  xend=lubridate::dmy(
                    c("31/12/2013", "1/5/2014", "31/12/2014", "1/5/2015", "31/12/2015")),
                  compartment="Susceptibles"),
                data.frame(
                  xbeg=lubridate::dmy(
                    c("13/11/2013", "1/3/2014", "1/9/2014", "1/3/2015", "1/9/2015")),
                  xend=lubridate::dmy(
                    c("31/12/2013", "1/5/2014", "31/12/2014", "1/5/2015", "31/12/2015")),
                  compartment="Infected"),
                data.frame(
                  xbeg=lubridate::dmy(
                    c("13/11/2013", "1/3/2014", "1/9/2014", "1/3/2015", "1/9/2015")),
                  xend=lubridate::dmy(
                    c("31/12/2013", "1/5/2014", "31/12/2014", "1/5/2015", "31/12/2015")),
                  compartment="Recovered")
              )

  int_data <- int_data %>%
                dplyr::left_join(
                  .,
                  date_iso_l %>%
                    dplyr::select(t, datesem),
                    by="t")

  # This chunk extracts the estimates and the credible intervals from calculations done in "other_functions.R"
  # The results of this chunk are already available in minimal data sets in /data (see README.md)

  # num_est <- lapply(1:length(estimate), function(x){
  #             return(
  #               estimate[[x]]$timeseries %>%
  #               dplyr::mutate(
  #                 sample=x))
  #           }) %>%
  #             do.call("rbind", .) %>%
  #           dplyr::group_by(sample) %>%
  #           dplyr::summarize(
  #             i=sum(i),
  #             i_v=sum(i_v),
  #             i_w=sum(i_w),
  #             i_w_v=sum(i_w_v)) %>%
  #           dplyr::mutate(
  #             red_v=i_v-i,
  #             red_w=i_w-i,
  #             red_w_v=i_w_v-i) %>%
  #           dplyr::ungroup() %>%
  #           dplyr::summarize(
  #             med_v=median(red_v),
  #             med_w=median(red_w),
  #             med_w_v=median(red_w_v),
  #             mean_v=mean(red_v),
  #             mean_w=mean(red_w),
  #             mean_w_v=mean(red_w_v),
  #             l_v=bayestestR::ci(red_v, method="HDI")$CI_low,
  #             l_w=bayestestR::ci(red_w, method="HDI")$CI_low,
  #             l_w_v=bayestestR::ci(red_w_v, method="HDI")$CI_low,
  #             up_v=bayestestR::ci(red_v, method="HDI")$CI_high,
  #             up_w=bayestestR::ci(red_w, method="HDI")$CI_high,
  #             up_w_v=bayestestR::ci(red_w_v, method="HDI")$CI_high,
  #             min_v=min(red_v),
  #             min_w=min(red_w),
  #             min_w_v=min(red_w_v),
  #             max_v=max(red_v),
  #             max_w=max(red_w),
  #             max_w_v=max(red_w_v)
  #             )

  load("PATH ON YOUR MACHINE/data/num_est.rda")

  box_d <- rbind(
          num_est %>%
            dplyr::select(
              min=min_v,
              med=med_v,
              mean=mean_v,
              max=max_v,
              l=l_v,
              up=up_v) %>%
            dplyr::mutate(
              type="No vaccination"),
          num_est %>%
            dplyr::select(
              min=min_w,
              med=med_w,
              mean=mean_w,
              max=max_w,
              l=l_w,
              up=up_w) %>%
            dplyr::mutate(
              type="No WASH"),
          num_est %>%
            dplyr::select(
              min=min_w_v,
              med=med_w_v,
              mean=mean_w_v,
              max=max_w_v,
              l=l_w_v,
              up=up_w_v) %>%
            dplyr::mutate(
              type="No vaccination\nand no WASH")
            ) %>%
          dplyr::mutate(
            type=factor(
              type,
              levels=c("No vaccination",
                  "No WASH",
                  "No vaccination\nand no WASH")))

  load("PATH ON YOUR MACHINE/data/vacci_res_c1.rda")
  vacci_opt <- vacci_res_c1 %>%
                dplyr::left_join(
                  .,
                  date_iso_l %>%
                    dplyr::select(t, datesem),
                    by="t")

  interv1_eps <- ggplot2::ggplot() +
        ggplot2::geom_rect(
          data=season_df,
          ggplot2::aes(
            xmin=xbeg,
            xmax=xend,
            ymin=-Inf,
            ymax=Inf),
          fill="grey",
          # alpha=0.2,
          color=NA) +
        ggplot2::geom_line(
          data=int_data,
          ggplot2::aes(
            x=datesem,
            y=value*100,
            color=type,
            size=type,
            group=type)) +
        ggplot2::scale_size_manual(
          name="",
          values=c(3,2.5,1.5,1),
          breaks=c("Intervention as\nit happened",
            "No vaccination",
            "No WASH",
            "No vaccination\nand no WASH")) +
        ggplot2::scale_color_manual(
          name="",
          values=viridis::viridis(n=4, option="D"),
          breaks=c("Intervention as\nit happened",
            "No vaccination",
            "No WASH",
            "No vaccination\nand no WASH")) +
        ggplot2::facet_wrap(.~compartment, ncol=1, scales="free_y") +
        ggplot2::geom_vline(
          xintercept=date_iso$datesem[1],
          color="darkviolet",
          linetype="dashed") +
        ggplot2::geom_vline(
          xintercept=date_iso$datesem[32],
          color="darkviolet",
          linetype="dashed") +
        ggplot2::geom_vline(
          xintercept=date_iso$datesem[35],
          color="darkviolet",
          linetype="dashed") +
        ggplot2::ylab("% of the population") +
        ggplot2::xlab("") +
        ggplot2::scale_x_date(
          labels=scales::label_date_short("%b-%Y"),
          breaks="years",
          limits=range(date_iso_l$datesem)) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          legend.position=c(0.15, 0.9),
          legend.margin = margin(t=-10),
          legend.spacing.x = unit(0, "mm"),
          legend.box.spacing = unit(0, "mm"),
          legend.text=ggplot2::element_text(size=15),
          axis.title.y=ggplot2::element_text(size = 26),
          axis.text.y=ggplot2::element_text(size = 26),
          axis.title.x=ggplot2::element_text(size = 26),
          axis.text.x=ggplot2::element_text(size = 26),
          strip.text = ggplot2::element_text(
            size = 20,
            margin = ggplot2::margin(t=12, b=12)))

  interv2 <- ggplot2::ggplot() +
        ggplot2::geom_violin(
          data=cum_df ,
          ggplot2::aes(y=value, x=type, color=type),
          fill=NA,
          size=2) +
        ggplot2::geom_errorbar(
          data=box_d,
          ggplot2::aes(ymin=l, ymax=up, x=type),
          width=0.1,
          size=1,
          color="black") +
        ggplot2::geom_point(
          data=box_d,
          ggplot2::aes(y=med, x=type, fill=type),
          shape=21,
          size=3,
          stroke=1,
          color="black") +
        ggplot2::geom_point(
          data=box_d,
          ggplot2::aes(y=mean, x=type, fill=type),
          shape="-",
          size=10,
          stroke=1,
          color="black") +
        ggplot2::scale_color_manual(
          name="",
          values=viridis::viridis(n=4, option="D")[-1],
          breaks=c(
            "No vaccination",
            "No WASH",
            "No vaccination\nand no WASH")) +
        ggplot2::scale_fill_manual(
          name="",
          values=viridis::viridis(n=4, option="D")[-1],
          breaks=c(
            "No vaccination",
            "No WASH",
            "No vaccination\nand no WASH")) +
        ggplot2::ylim(0, NA) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          legend.position="none",
          axis.title.y=ggplot2::element_text(size = 26),
          axis.text.y=ggplot2::element_text(size = 26),
          axis.title.x=ggplot2::element_text(size = 26),
          axis.text.x=ggplot2::element_text(size = 26),
          strip.text = ggplot2::element_text(size = 26)) +
        ggplot2::ylab("Number of additionnal cases compared\nto the original intervention") +
        ggplot2::xlab("")

  vacci_time <- ggplot2::ggplot(
          data=vacci_opt %>%
            dplyr::mutate(
              col=ifelse(l<=0 & u>=0, NA, "black"))) +
          ggplot2::geom_tile(
            ggplot2::aes(x=datesem, y=100*N/262963, fill=avoid),
            color="black") +
        ggplot2::geom_vline(
          xintercept=date_iso$datesem[1],
          color="darkviolet",
          linetype="dashed") +
        ggplot2::geom_vline(
          xintercept=date_iso$datesem[32],
          color="darkviolet",
          linetype="dashed") +
        ggplot2::geom_vline(
          xintercept=date_iso$datesem[35],
          color="darkviolet",
          linetype="dashed") +
          ggplot2::scale_fill_gradient2(
            name="Cases avoided",
            high="gray30",
            limits=c(0, 10000)) +
        ggplot2::scale_x_date(
          labels=scales::label_date_short("%b-%Y"),
          breaks="years") +
        ggplot2::xlab("") +
        ggplot2::ylab("% of vaccinated population") +
        ggplot2::theme_bw() +
        ggplot2::theme(
          legend.position="top",
          legend.text=ggplot2::element_text(size=15),
          axis.title.y = ggplot2::element_text(size = 26),
          axis.text.y = ggplot2::element_text(size = 26),
          axis.title.x = ggplot2::element_text(size = 26),
          axis.text.x = ggplot2::element_text(size = 26))

  library(patchwork)

  ggplot2::ggsave(
    "PATH TO YOUR FOLDER/fig3.eps",
    interv1_eps +
      (interv2 +
        vacci_time +
        plot_layout(nrow=2)) +
      plot_layout(
        widths=c(3,2)),
    width=16,
    height=10)
  # The final figure available in the manuscript was then improved in adobe illustrator

}

fig4 <- function(){
  # This chunk extracts the estimates and the credible intervals from calculations done in "other_functions.R"
  # The results of this chunk are already available in minimal data sets in /data (see README.md)

  # No env and no contamination
  # prop_env <- lapply(1:length(contrib), function(x){
  #   return(
  #     contrib[[x]]$prop %>%
  #     dplyr::filter(t %in% 0:118) %>%
  #     dplyr::mutate(
  #       sample=x))
  # }) %>%
  #   do.call("rbind", .)

  date_iso <- data.frame(
                year=c(rep(2013, length(47:52)), rep(2014, 52), rep(2015, 53), rep(2016, 7)),
                week=c(47:52, 1:52, 1:53, 1:7)) %>%
                dplyr::mutate(
                  datesem=ISOweek::ISOweek2date(paste0(year, "-W", stringr::str_pad(week, 2, side="left", pad = "0"), "-3")))

  # foi_env <- prop_env %>%
  #             dplyr::group_by(t) %>%
  #             dplyr::summarize(
  #               l_h=bayestestR::ci(foi_h, method="HDI")$CI_low,
  #               l_e=bayestestR::ci(foi_e, method="HDI")$CI_low,
  #               up_h=bayestestR::ci(foi_h, method="HDI")$CI_high,
  #               up_e=bayestestR::ci(foi_e, method="HDI")$CI_high,
  #               foi_h=mean(foi_h),
  #               foi_e=mean(foi_e)) %>%
  #             dplyr::ungroup() %>%
  #             dplyr::left_join(
  #               .,
  #               date_iso_l %>%
  #                 dplyr::select(t, datesem),
  #               by="t")

  # growth_env  <- prop_env %>%
  #             dplyr::group_by(t) %>%
  #             dplyr::summarize(
  #               l=bayestestR::ci(-net_decay, method="HDI")$CI_low,
  #               up=bayestestR::ci(-net_decay, method="HDI")$CI_high,
  #               growth=-mean(net_decay)) %>%
  #             dplyr::ungroup() %>%
  #             dplyr::left_join(
  #               .,
  #               date_iso_l %>%
  #                 dplyr::select(t, datesem),
  #               by="t")

  # env_data <- rbind(
  #       prop_env %>%
  #         dplyr::select(
  #           t, Susceptibles=s, Infected=i) %>%
  #         dplyr::mutate(Recovered=1-Susceptibles-Infected) %>%
  #         tidyr::gather(., compartment, value, Susceptibles:Recovered) %>%
  #         dplyr::group_by(t, compartment) %>%
  #         dplyr::summarize(
  #           value=mean(value)) %>%
  #         dplyr::ungroup() %>%
  #         dplyr::mutate(type="No vacc.\nand no WASH"),
  #       prop_env %>%
  #         dplyr::select(
  #           t, Susceptibles=s_h, Infected=i_h) %>%
  #         dplyr::mutate(Recovered=1-Susceptibles-Infected) %>%
  #         tidyr::gather(., compartment, value, Susceptibles:Recovered) %>%
  #         dplyr::group_by(t, compartment) %>%
  #         dplyr::summarize(
  #           value=mean(value)) %>%
  #         dplyr::ungroup() %>%
  #         dplyr::mutate(type="No seas.\nmigr."),
  #       prop_env %>%
  #         dplyr::select(
  #           t, Susceptibles=s_e, Infected=i_e) %>%
  #         dplyr::mutate(Recovered=1-Susceptibles-Infected) %>%
  #         tidyr::gather(., compartment, value, Susceptibles:Recovered) %>%
  #         dplyr::group_by(t, compartment) %>%
  #         dplyr::summarize(
  #           value=mean(value)) %>%
  #         dplyr::ungroup() %>%
  #         dplyr::mutate(type="No env.\nexp."),
  #       prop_env %>%
  #         dplyr::select(
  #           t, Susceptibles=s_c, Infected=i_c) %>%
  #         dplyr::mutate(Recovered=1-Susceptibles-Infected) %>%
  #         tidyr::gather(., compartment, value, Susceptibles:Recovered) %>%
  #         dplyr::group_by(t, compartment) %>%
  #         dplyr::summarize(
  #           value=mean(value)) %>%
  #         dplyr::ungroup() %>%
  #         dplyr::mutate(type="No env.\ncont.")
  #       ) %>%
  #       dplyr::mutate(
  #         compartment=factor(
  #           compartment,
  #           levels=c("Infected", "Recovered", "Susceptibles")),
  #         type=factor(
  #           type,
  #           levels=c("No vacc.\nand no WASH",
  #                 "No seas.\nmigr.",
  #                 "No env.\nexp.",
  #                 "No env.\ncont."))) %>%
  #       dplyr::left_join(
  #         .,
  #         date_iso_l %>%
  #           dplyr::select(t, datesem),
  #           by="t")

  load("PATH ON YOUR MACHINE/data/foi_env.rda")
  load("PATH ON YOUR MACHINE/data/growth_env.rda")
  load("PATH ON YOUR MACHINE/data/env_data.rda")

  date_iso_l <- rbind(
                  data.frame(
                    year=2013,
                    week=46,
                    datesem=lubridate::dmy("13/11/2013")),
                  date_iso) %>%
                dplyr::mutate(
                  t=0:118)

  season_df <- rbind(
                data.frame(
                  xbeg=lubridate::dmy(
                    c("13/11/2013", "1/3/2014", "1/9/2014", "1/3/2015", "1/9/2015")),
                  xend=lubridate::dmy(
                    c("31/12/2013", "1/5/2014", "31/12/2014", "1/5/2015", "31/12/2015")),
                  compartment="Susceptibles"),
                data.frame(
                  xbeg=lubridate::dmy(
                    c("13/11/2013", "1/3/2014", "1/9/2014", "1/3/2015", "1/9/2015")),
                  xend=lubridate::dmy(
                    c("31/12/2013", "1/5/2014", "31/12/2014", "1/5/2015", "31/12/2015")),
                  compartment="Infected"),
                data.frame(
                  xbeg=lubridate::dmy(
                    c("13/11/2013", "1/3/2014", "1/9/2014", "1/3/2015", "1/9/2015")),
                  xend=lubridate::dmy(
                    c("31/12/2013", "1/5/2014", "31/12/2014", "1/5/2015", "31/12/2015")),
                  compartment="Recovered")
              )

  foi_env <- rbind(
              foi_env %>%
                dplyr::select(datesem, foi=foi_e, l=l_e, up=up_e) %>%
                dplyr::mutate(type="Environment"),
              foi_env %>%
                dplyr::select(datesem, foi=foi_h, l=l_h, up=up_h) %>%
                dplyr::mutate(type="Interhuman")
              )

  # This chunk extracts the estimates and the credible intervals from calculations done in "other_functions.R"
  # The results of this chunk are already available in minimal data sets in /data (see README.md)

  # cum_df2 <- lapply(1:length(contrib), function(x){
  #           return(
  #             contrib[[x]]$timeseries %>%
  #             dplyr::mutate(
  #               sample=x))
  #         }) %>%
  #         do.call("rbind", .) %>%
  #         dplyr::group_by(sample) %>%
  #         dplyr::summarize(
  #             red_h=sum(h),
  #             red_e=sum(e),
  #             red_c=sum(c)) %>%
  #         dplyr::ungroup() %>%
  #         tidyr::gather(., type, value, red_h:red_c) %>%
  #         dplyr::mutate(
  #           type=factor(
  #             type,
  #             levels=c("red_h", "red_e", "red_c"),
  #             labels=c("No seas.\nmigr.",
  #                 "No env.\nexp.",
  #                 "No env.\ncont.")))

  load("PATH ON YOUR MACHINE/data/cum_df2.rda")

  # This chunk extracts the estimates and the credible intervals from calculations done in "other_functions.R"
  # The results of this chunk are already available in minimal data sets in /data (see README.md)

  # num_env <- lapply(1:length(contrib), function(x){
  #               return(
  #                 contrib[[x]]$timeseries %>%
  #                 dplyr::select(t, h, e, c) %>%
  #                 dplyr::mutate(
  #                   sample=x))
  #             }) %>%
  #           do.call("rbind", .) %>%
  #           dplyr::group_by(sample) %>%
  #           dplyr::summarize(
  #             red_h=sum(h),
  #             red_e=sum(e),
  #             red_c=sum(c)) %>%
  #           dplyr::ungroup() %>%
  #           dplyr::summarize(
  #             med_h=median(red_h),
  #             med_e=median(red_e),
  #             med_c=median(red_c),
  #             mean_h=mean(red_h),
  #             mean_e=mean(red_e),
  #             mean_c=mean(red_c),
  #             l_h=bayestestR::ci(red_h, method="HDI")$CI_low,
  #             l_e=bayestestR::ci(red_e, method="HDI")$CI_low,
  #             l_c=bayestestR::ci(red_c, method="HDI")$CI_low,
  #             up_h=bayestestR::ci(red_h, method="HDI")$CI_high,
  #             up_e=bayestestR::ci(red_e, method="HDI")$CI_high,
  #             up_c=bayestestR::ci(red_c, method="HDI")$CI_high,
  #             min_h=min(red_h),
  #             min_e=min(red_e),
  #             min_c=min(red_c),
  #             max_h=max(red_h),
  #             max_e=max(red_e),
  #             max_c=max(red_c)
  #             )

  load("PATH ON YOUR MACHINE/data/num_env.rda")

  box_d <- rbind(
          num_env %>%
            dplyr::select(
              min=min_e,
              med=med_e,
              mean=mean_e,
              max=max_e,
              l=l_e,
              up=up_e) %>%
            dplyr::mutate(
              type="No env.\nexp."),
          num_env %>%
            dplyr::select(
              min=min_c,
              med=med_c,
              mean=mean_c,
              max=max_c,
              l=l_c,
              up=up_c) %>%
            dplyr::mutate(
              type="No env.\ncont.")
            ) %>%
          dplyr::mutate(
            type=factor(
              type,
              levels=c(
                  "No env.\nexp.",
                  "No env.\ncont.")))

  # We extract the average values for beta_h and beta_e
  # We directly provide them below

  # env_sample <- build_sample2(name=paste0(c("theta1_", "theta2_", "theta3_"), "0_8_s3_migrN"))

  # theta1_sample <- env_sample$theta1 %>%
  #                 setNames(., c("alpha0", "alpha2",
  #                   "alpha4", "alpha5", "beta_env",
  #                   "beta_wash", "epsilon", "cont", "decay",
  #                   "lambda", "lambda_c"))

  # beta_h <- exp(mean(theta1_sample$alpha0))
  # beta_e <- mean(theta1_sample$beta_env)

  beta_h <- 1.361667
  beta_e <- 0.8437898

  human_env_eps <- ggplot2::ggplot() +
              ggplot2::geom_rect(
                data=season_df,
                ggplot2::aes(
                  xmin=xbeg,
                  xmax=xend,
                  ymin=-Inf,
                  ymax=Inf),
                fill="grey",
                # alpha=0.2,
                color=NA) +
              ggplot2::geom_line(
                data=env_data,
                ggplot2::aes(
                  x=datesem,
                  y=value*100,
                  color=type,
                  size=type,
                  group=type)) +
              ggplot2::scale_size_manual(
                name="",
                values=rev(c(1,1.5,2,3)),
                breaks=c("No vacc.\nand no WASH",
                  "No seas.\nmigr.",
                  "No env.\nexp.",
                  "No env.\ncont.")) +
              ggplot2::scale_color_manual(
                name="",
                values=rcartocolor::carto_pal(4, "Safe"),
                breaks=c("No vacc.\nand no WASH",
                  "No seas.\nmigr.",
                  "No env.\nexp.",
                  "No env.\ncont.")) +
              ggplot2::facet_wrap(.~compartment, ncol=1, scales="free_y") +
              ggplot2::guides(
                col=ggplot2::guide_legend(ncol=2)) +
              ggplot2::ylab("% of the population") +
              ggplot2::xlab("") +
              ggplot2::scale_x_date(
                labels=scales::label_date_short("%b-%Y"),
                breaks="years",
                limits=range(date_iso_l$datesem)) +
              ggplot2::theme_bw() +
              ggplot2::theme(
                # legend.position="bottom",
                legend.position=c(0.15, 0.9),
                legend.margin = margin(t=-10),
                # legend.spacing.y = unit(0, "mm"),
                legend.spacing.x = unit(0, "mm"),
                legend.box.spacing = unit(0, "mm"),
                legend.text=ggplot2::element_text(size=15),
                axis.title.y=ggplot2::element_text(size = 26),
                axis.text.y=ggplot2::element_text(size = 26),
                axis.title.x=ggplot2::element_text(size = 26),
                axis.text.x=ggplot2::element_text(size = 26),
                strip.text = ggplot2::element_text(
                  size = 20,
                  margin = ggplot2::margin(t=12, b=12)))

  y_upper_limit <- diff(range(cum_df2$value)) * 0.05 + max(cum_df2$value)
  y_lower_limit <- 0 - diff(range(cum_df2$value)) * 0.05

  human_env2 <- ggplot2::ggplot() +
        ggplot2::geom_violin(
          data=cum_df2,
          ggplot2::aes(y=value, x=type, color=type),
          fill=NA,
          size=2) +
        ggplot2::geom_errorbar(
          data=box_d,
          ggplot2::aes(ymin=l, ymax=up, x=type),
          width=0.1,
          size=1,
          color="black") +
        ggplot2::geom_point(
          data=box_d,
          ggplot2::aes(y=med, x=type, fill=type),
          shape=21,
          size=3,
          stroke=1,
          color="black") +
        ggplot2::geom_point(
          data=box_d,
          ggplot2::aes(y=mean, x=type, fill=type),
          shape="-",
          size=10,
          stroke=1,
          color="black") +
        ggplot2::scale_color_manual(
          name="",
          values=rcartocolor::carto_pal(4, "Safe")[-1],
          breaks=c(
            "No seas.\nmigr.",
            "No env.\nexp.",
            "No env.\ncont.")) +
        ggplot2::scale_fill_manual(
          name="",
          values=rcartocolor::carto_pal(4, "Safe")[-1],
          breaks=c(
            "No seas.\nmigr.",
            "No env.\nexp.",
            "No env.\ncont.")) +
        ggplot2::facet_grid(
          .~type,
          # nrow=1,
          scales="free_x",
          space="free_x") +
        ggplot2::geom_hline(yintercept = y_upper_limit) +
        ggplot2::scale_y_continuous(
          labels = ~ format(.x, scientific = FALSE),
          sec.axis = ggplot2::dup_axis(breaks = NULL),
          expand = c(0, 0)) +
        ggplot2::coord_cartesian(ylim = c(y_lower_limit, y_upper_limit)) +
        ggplot2::theme_classic() +
        # ggplot2::theme_bw() +
        ggplot2::theme(
          legend.position="none",
          panel.background = ggplot2::element_rect(fill = "white", 
            colour = NA),
          panel.grid = ggplot2::element_line(colour = "grey92"), 
          panel.grid.minor = ggplot2::element_line(linewidth = rel(0.5)),
          panel.grid.major = ggplot2::element_line(linewidth = rel(0.75)),
          axis.title.y.right = ggplot2::element_blank(),
          axis.text.y.right = ggplot2::element_blank(),
          # axis.ticks.y = ggplot2::element_blank(),
          axis.title.y=ggplot2::element_text(size = 26),
          axis.text.y=ggplot2::element_text(size = 26),
          axis.title.x=ggplot2::element_text(size = 26),
          axis.text.x=ggplot2::element_text(size = 26),
          strip.text = ggplot2::element_blank(),
          # panel.border = ggplot2::element_blank(),
          panel.margin = unit(0, "mm"),
          strip.background = ggplot2::element_rect(linewidth = 0.5)) +
        ggplot2::ylab("Avoided cases") +
        ggplot2::xlab("")

  foi_fig_eps <- ggplot2::ggplot() +
              ggplot2::geom_rect(
                data=season_df,
                ggplot2::aes(
                  xmin=xbeg,
                  xmax=xend,
                  ymin=-Inf,
                  ymax=Inf),
                fill="grey",
                color=NA) +
              ggplot2::geom_ribbon(
                data=foi_env,
                ggplot2::aes(
                  x=datesem,
                  ymin=l,
                  ymax=up,
                  fill=type),
                color=NA) +
              ggplot2::geom_line(
                data=foi_env,
                ggplot2::aes(
                  x=datesem,
                  y=foi,
                  color=type),
                size=1) +
              ggplot2::geom_hline(
                data=data.frame(
                  value=c(beta_h, beta_e),
                  type=c("Interhuman", "Environment")),
                ggplot2::aes(
                  yintercept=value,
                  color=type),
                size=1,
                linetype="dashed") +
              ggplot2::scale_color_brewer(
                name="",
                palette="Set1",
                labels=c("Environmental", "Interhuman")) +
              ggplot2::scale_fill_brewer(
                name="",
                palette="Set1",
                labels=c("Environmental", "Interhuman")) +
              ggplot2::ylab("FOI") +
              ggplot2::xlab("") +
              ggplot2::scale_x_date(
                labels=scales::label_date_short("%b-%Y"),
                breaks="years",
                limits=range(date_iso_l$datesem)) +
              ggplot2::theme_bw() +
              ggplot2::theme(
                legend.position=c(0.2, 0.8),
                legend.margin = margin(t=-10),
                legend.spacing.y = unit(0, "mm"),
                legend.spacing.x = unit(0, "mm"),
                legend.box.spacing = unit(0, "mm"),
                legend.text=ggplot2::element_text(size = 20),
                axis.title.y=ggplot2::element_text(size = 26),
                axis.text.y=ggplot2::element_text(size = 26),
                axis.title.x=ggplot2::element_blank(),
                axis.text.x=ggplot2::element_blank(),
                axis.ticks.x=ggplot2::element_blank(),
                plot.margin=ggplot2::margin(r=0, b=0))

  growth_fig_eps <- ggplot2::ggplot() +
                ggplot2::geom_ribbon(
                  data=growth_env,
                  ggplot2::aes(
                    x=datesem,
                    ymin=l,
                    ymax=up),
                  fill="gray25",
                  # alpha=0.2,
                  color=NA) +
                ggplot2::geom_line(
                  data=growth_env,
                  ggplot2::aes(
                    x=datesem,
                    y=growth),
                  color="gray35",
                  size=1) +
                ggplot2::geom_hline(
                  yintercept=0,
                  linetype="dashed",
                  color="black",
                  size=1) +
                ggplot2::ylim(-0.7, 1.6) +
                ggplot2::ylab("Growth rate") +
                ggplot2::xlab("") +
                ggplot2::scale_x_date(
                  labels=scales::label_date_short("%b-%Y"),
                  breaks="years",
                  limits=range(date_iso_l$datesem)) +
                ggplot2::theme_bw() +
                ggplot2::theme(
                  axis.title.y=ggplot2::element_text(size = 26),
                  axis.text.y=ggplot2::element_text(size = 26),
                  axis.title.x=ggplot2::element_text(size = 26),
                  axis.text.x=ggplot2::element_text(size = 26),
                  strip.text = ggplot2::element_text(size = 26),
                  plot.margin=ggplot2::margin(r=0, t=0))

  growth_dist <- ggplot2::ggplot() +
                ggplot2::geom_violin(
                  data=growth_env,
                  ggplot2::aes(
                    x=1,
                    y=growth),
                  fill=NA,
                  color="darkgrey",
                  width=1,
                  # size=2,
                  linewidth=2) +
                ggplot2::geom_boxplot(
                  data=growth_env,
                  ggplot2::aes(
                    x=1,
                    y=growth),
                  fill="darkgrey",
                  size=0.2,
                  width=0.2) +
                ggplot2::geom_hline(
                  yintercept=0,
                  linetype="dashed",
                  color="black",
                  size=1) +
                ggplot2::ylim(-0.7, 1.6) +
                ggplot2::xlim(0.45, 1.55) +
                ggplot2::ylab("") +
                ggplot2::xlab("") +
                ggplot2::theme_bw() +
                ggplot2::theme(
                  axis.title.y=ggplot2::element_blank(),
                  axis.text.y=ggplot2::element_blank(),
                  axis.ticks.y=ggplot2::element_blank(),
                  axis.title.x=ggplot2::element_blank(),
                  axis.text.x=ggplot2::element_blank(),
                  axis.ticks.x=ggplot2::element_blank(),
                  plot.margin=ggplot2::margin(l=0, t=0))

  library(patchwork)

  ggplot2::ggsave(
    "PATH TO YOUR FOLDER/fig4.eps",
    human_env_eps +
      (human_env2 +
        (foi_fig_eps +
          plot_spacer() +
          plot_layout(ncol=2, widths= c(7, 1.5))) +
        (growth_fig_eps +
          growth_dist +
          plot_layout(ncol=2, widths= c(7, 1.5))) +
        plot_layout(nrow=3, heights=c(1, 0.5, 0.5))) +
      plot_layout(
        widths=c(3,2)),
    width=16,
    height=10)
  # The final figure available in the manuscript was then improved in adobe illustrator

}

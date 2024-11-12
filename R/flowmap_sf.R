#' Create sf for flowmap
#'
#' @param flowdat Input dataframe. See details below.
#' @param od As an alternative to \code{flowdat}, dataframe with the origin-destination pairs and the flow between them.  Must contain the columns o, d, value. \code{nodes} must be provided as well. See details below.
#' @param nodes As an alternative to \code{flowdat}, a dataframe with the nodes of the network. Must contain the columns name, x, y. See details below.
#' @param crs A spatial reference system (CRS) code of tha coordinates passed in \code{flowdat} or \code{od}.
#' @param k_nodes Number of clusters to group nodes into. If defined, nodes will be clustered hierarchically based on spatial proximity. By default, no clustering will be applied.
#' @param node_buffer_factor Controls the distance between the nodes and the edges ( in multiple of the nodes' radii).
#' @param node_radius_factor Controls the size of the nodes.
#' @param node_fill_factor Controls the downscaling of the fill of the nodes ( as to not outshine the edges ).
#' @param edge_offset_factor Controls the distance between the parallel arrows.
#' @param edge_width_factor Controls the width of the edges.
#' @param outline_linewidth The linewidth of the outline of the arrows.
#' @param alpha Opacity of the edges.
#' @param nodes_alpha Opacity of the nodes.
#' @param outline_col Color of the outline of the edges.
#' @param arrow_point_angle Controls the pointiness of the edges.
#'
#' @details
#' The function requires as inputs a dataframe \code{flowdat} which contains for every combination of two nodes a and b the coordinates of these nodes as well as the intensity of flow between those nodes in both directions (a to b, b to a). The dataframe should have the following columns:
#' \itemize{
#'  \item \strong{id_a:} The unique id of node a
#'  \item \strong{id_b:} The unique id of node b
#'  \item \strong{xa:} The x coordinate of node a
#'  \item \strong{ya:} The y coordinate of node a
#'  \item \strong{xb:} The x coordinate of node b
#'  \item \strong{yb:} The y coordinate of node b
#'  \item \strong{flow_ab:} The intensity of flow from node a to node b
#'  \item \strong{flow_ba:} The intensity of flow from node b to node a
#' }
#'
#' Alternatively, the function can take as input a dataframe \code{od} which contains the origin-destination pairs and the flow between them. The dataframe should have the following columns:
#' \itemize{
#' \item \strong{o:} The unique id of the origin node
#' \item \strong{d:} The unique id of the destination node
#' \item \strong{value:} The intensity of flow between the origin and destination
#' }
#' In this case, the function also requires a dataframe \code{nodes} which contains the coordinates of the nodes. The dataframe should have the following columns:
#' \itemize{
#' \item \strong{name:} The unique id of the node
#' \item \strong{x:} The x coordinate of the node
#' \item \strong{y:} The y coordinate of the node
#' }
#'
#' The function will impose coord_equal() on the ggplot.
#'
#' Inspired by \href{https://flowmap.gl/}{flowmap.gl}.
#'
#' @importFrom dplyr mutate select left_join summarize group_by ungroup bind_rows n arrange group_split n_distinct across where filter_all coalesce all_vars
#' @importFrom tidyr pivot_longer separate
#' @importFrom forcats fct_reorder
#' @importFrom ggplot2 ggplot geom_polygon annotate
#' @return The ggplot with an additional polygon layer for the flow arrows and an additional polygon layer for the nodes
#' @export
#'
#' @examples
#' data <- flowmapper::CH_migration_data
#' cantons <- flowmapper::cantons
#' sf::st_crs(cantons) <- 3857
#' flows_nodes <- flowmap_sf(flowdat = data, nodes = cantons, crs = 3857)
flowmap_sf <- function(flowdat=NULL,od=NULL,nodes=NULL,crs=NULL,outline_linewidth=0.01,alpha=0.8,nodes_alpha=0.8,outline_col="black",k_nodes=NULL,node_buffer_factor = 1.2, node_radius_factor = 1, edge_offset_factor = 1, node_fill_factor = NULL, edge_width_factor = 1.2, arrow_point_angle = 45){

  if(!is.null(od) & is.null(nodes)){
    stop("When providing data in od format, nodes (a dataframe with the nodes coordiantes) must be defined as well.")
  }

  if(!is.null(od) & is.null(nodes)){
    stop("When providing nodes, od must be defined as well.")
  }

  if(!is.null(flowdat) & !is.null(od)){
    stop("Only one of flowdat or od can be defined.")
  }

  if(is.null(flowdat) & is.null(od)){
    stop("Either flowdat or od must be defined.")
  }

  if(is.null(flowdat) & !is.null(od) & !is.null(nodes)){
    #force convert columns o and d to characters
    od <- od |> mutate(o=as.character(o),d=as.character(d))
    # force convert column name to character
    nodes <- nodes |> mutate(name=as.character(name))
    flowdat <- util_data_flow_to_flowdat(nodes,od)
  }

  if(inherits(flowdat,"grouped_df")){
    flowdat <- flowdat |> dplyr::ungroup()
  }
  # Some checks
  if(arrow_point_angle>75){arrow_point_angle <- 75; warning("Warning. arrow_point_angle cannot exceed 75 degrees")}
  if(arrow_point_angle<20){arrow_point_angle <- 20; warning("Warning. arrow_point_angle cannot be lower than 20 degrees")}

  if(!is.null(k_nodes)){
    flowdat <- hca_flowdat(flowdat,k_nodes)
  }else{
    if(n_distinct(flowdat$id_a)>50){
      warning("Number of flows very high. Consider setting k_nodes to cluster nodes.")
    }
  }

  flowdat <- flowdat |> dplyr::filter(id_a!=id_b)

  #force convert factor columns to character
  flowdat <- flowdat |> dplyr::mutate(dplyr::across(dplyr::where(is.factor), as.character))

  # remove any row with a missing value in any column
  flowdat <- flowdat |> dplyr::filter_all(dplyr::all_vars(!is.na(.)))

  nodes <-
    bind_rows(
      flowdat |>
        mutate(flow=flow_ab+flow_ba) |>
        select(id=id_a,x=xa,y=ya,flow),
      flowdat |>
        mutate(flow=flow_ab+flow_ba) |>
        select(id=id_b,x=xb,y=yb,flow)
    ) |>
    group_by(id,x,y) |>
    summarize(
      flowsum=sum(flow),
      radius=sqrt(flowsum/pi),
      n_edges=n(),
    ) |>
    ungroup()


  if(is.null(node_fill_factor)){
    node_fill_factor <- max(c(flowdat$flow_ab,flowdat$flow_ba))/max(nodes$flowsum)
  }

  xrange=max(nodes$x)-min(nodes$x)
  yrange=max(nodes$y)-min(nodes$y)
  # width of the widest arrow
  maxwidth <- 0.05*max(xrange,yrange)*edge_width_factor
  # circle size
  max_circle_size <- 0.02*max(xrange,yrange)*node_radius_factor
  # value by which the two arrows are separated
  off <- 0.005*max(xrange,yrange)*edge_offset_factor

  # create group id
  flowdat <-
    flowdat |>
    mutate(group=paste0(id_a, " - ", id_b))

  # calculate the adjusted radius for the nodes
  nodes <-
    nodes |>
    mutate(
      adj_radius=max_circle_size*radius/max(radius))

  # compute global scaling
  max_flow <- max(max(flowdat$flow_ab),max(flowdat$flow_ba))

  ###########
  flowdat_ab <-
    flowdat |>
    # add the information about the radius of the nodes
    left_join(
      nodes |> select(id,adj_radius_b=adj_radius),
      by=c("id_a"="id")
    )|>
    left_join(
      nodes |> select(id,adj_radius_a=adj_radius),
      by=c("id_b"="id")
    ) |>
    mutate(width=maxwidth*flow_ab/max_flow,
           e = atan2(xb-xa,yb-ya),
           #apply shortening
           xa = xa+sin(e)*adj_radius_b*node_buffer_factor,
           ya = ya+cos(e)*adj_radius_b*node_buffer_factor,
           xb = xb-sin(e)*adj_radius_a*node_buffer_factor,
           yb = yb-cos(e)*adj_radius_a*node_buffer_factor,

           a_b = sqrt((xb-xa)^2+(yb-ya)^2),
           m_c = width*tan((90-arrow_point_angle)*(pi/180)),
           a_m = a_b - m_c,
           xm = xa+(xb-xa)*(a_m/a_b), # coordinates of B
           ym = ya+(yb-ya)*(a_m/a_b), # coordinates of B

           xe=cos(e)*width+xm,
           ye=sin(e)*-width+ym,
           xf=cos(e)*width/2+xm,
           yf=sin(e)*-width/2+ym,
           xg=cos(e)*width/2+xa,
           yg=sin(e)*-width/2+ya,
           # calculate offsets
           off_x=cos(e)*off/2,
           off_y=sin(e)*-off/2,
           # apply offsets
           xa=xa+off_x,xb=xb+off_x,xe=xe+off_x,xm=xm+off_x,xf=xf+off_x,xg=xg+off_x,
           ya=ya+off_y,yb=yb+off_y,ye=ye+off_y,ym=ym+off_y,yf=yf+off_y,yg=yg+off_y,
    )

  flowdat_ba <-
    flowdat  |>
    # add the information about the radius of the nodes
    left_join(
      nodes |> select(id,adj_radius_b=adj_radius),
      by=c("id_a"="id")
    )|>
    left_join(
      nodes |> select(id,adj_radius_a=adj_radius),
      by=c("id_b"="id")
    ) |>
    mutate(width=maxwidth*flow_ba/max_flow,
           # swap o and d coords, keep all else the same
           xz=xa,xa=xb,xb=xz,xz=NULL,
           yz=ya,ya=yb,yb=yz,yz=NULL,
           # end swap

           e = atan2(xb-xa,yb-ya),

           #apply shortening
           xa = xa+sin(e)*adj_radius_a*node_buffer_factor,
           ya = ya+cos(e)*adj_radius_a*node_buffer_factor,
           xb = xb-sin(e)*adj_radius_b*node_buffer_factor,
           yb = yb-cos(e)*adj_radius_b*node_buffer_factor,


           a_b = sqrt((xb-xa)^2+(yb-ya)^2),
           m_c = width*tan((90-arrow_point_angle)*(pi/180)),
           a_m = a_b - m_c,
           xm = xa+(xb-xa)*(a_m/a_b), # coordinates of m
           ym = ya+(yb-ya)*(a_m/a_b), # coordinates of m

           xe=cos(e)*width+xm,
           ye=sin(e)*-width+ym,
           xf=cos(e)*width/2+xm,
           yf=sin(e)*-width/2+ym,
           xg=cos(e)*width/2+xa,
           yg=sin(e)*-width/2+ya,

           # calculate offsets
           off_x=cos(e)*off/2,
           off_y=sin(e)*-off/2,
           # apply offsets
           xa=xa+off_x,xb=xb+off_x,xe=xe+off_x,xm=xm+off_x,xf=xf+off_x,xg=xg+off_x,
           ya=ya+off_y,yb=yb+off_y,ye=ye+off_y,ym=ym+off_y,yf=yf+off_y,yg=yg+off_y,
    )

  #calculate sum of flows for each route (for ordering the arrows drawing order)
  flowsums <-
    flowdat_ab |>
    select(group,flow_ab) |>
    left_join(flowdat_ba|>
                select(group,flow_ba),by="group") |>
    mutate(flowsum=flow_ab+flow_ba) |>
    select(group,flowsum)

  flows <-
    bind_rows(
      flowdat |>
        mutate(group = paste0(group,"_ab")) |>
        select(group,flow=flow_ab),
      flowdat |>
        mutate(group = paste0(group,"_ba")) |>
        select(group,flow=flow_ba)
    )

  flows <-
    bind_rows(
      flowdat |>
        mutate(group = group) |>
        select(group,flow=flow_ab),
      flowdat |>
        separate(group,into=c("aux1","aux2"),sep = " - ") |>
        mutate(group = paste0(aux2," - ",aux1)) |>
        select(group,flow=flow_ba)
    )

  plot_df <-
    bind_rows(
      left_join(
        flowdat_ab |>
          select(group,ya,ym,yb,ye,yf,yg) |>
          tidyr::pivot_longer(-group) |>
          mutate(point=sub('.', '', name)) |>
          select(group,point,y=value),
        flowdat_ab |>
          select(group,xa,xm,xb,xe,xf,xg) |>
          tidyr::pivot_longer(-group) |>
          mutate(point=sub('.', '', name)) |>
          select(group,point,x=value),
        by=c("group","point")
      ) |>
        mutate(dir="ab"),
      left_join(
        flowdat_ba |>
          select(group,ya,ym,yb,ye,yf,yg) |>
          tidyr::pivot_longer(-group) |>
          mutate(point=sub('.', '', name)) |>
          select(group,point,y=value),
        flowdat_ba |>
          select(group,xa,xm,xb,xe,xf,xg) |>
          tidyr::pivot_longer(-group) |>
          mutate(point=sub('.', '', name)) |>
          select(group,point,x=value),
        by=c("group","point")
      ) |>
        mutate(dir="ba")
    ) |>
    # order by sum of flows to have bigger flows on top
    left_join(flowsums,by="group") |>
    arrange(flowsum) |>
    #
    # mutate(group=paste0(group,"_",dir))|>
    # create unique group name for each direction
    separate(group,into=c("aux1","aux2"),sep = " - ") |>
    mutate(group=ifelse(
      dir=="ab",
      paste0(aux1," - ", aux2),
      paste0(aux2," - ", aux1)
    ))|>
    # add the flows for each arrow (for fill)
    left_join(flows,by="group") |>
    mutate(
      label=paste0("Route: ",group,"\nFlow: ",flow),
      group=forcats::fct_reorder(group,flowsum)
    )

  plot_df <-
    plot_df |> group_by(label) |> filter(sum(flow)>0) |> ungroup()

  # for the nodes, calculate the flow to be visualized (as fill) from the node_fill_factor and their flow sum
  nodes <-
    nodes |>
    mutate(flow=flowsum*node_fill_factor)

  nodes_poly <-
    nodes |>
    group_by(id) |>
    group_split() |>
    lapply(function(x){
      coords <- get_circle_coords(c(x$x,x$y),x$adj_radius)
      coords$group=x$id
      coords$flowsum=x$flowsum
      coords$flow=x$flow
      return(coords)
    }) |>
    bind_rows()|>
    mutate(label = paste0("Node: ",group,"\nTotal Flow: ",flowsum))

  # convert flows to sf
  flows_sf <- plot_df |> 
    arrange(group, point) |> 
    group_by(group) |> 
    summarize(
      flowsum = first(flowsum),
      flow = first(flow),
      label = first(label),
      geometry = sf::st_sfc(
        sf::st_polygon(list(rbind(cbind(x, y), cbind(x, y)[1, ]))))
    ) |> 
    sf::st_as_sf(crs = crs) # Set CRS as appropriate for your data

  
  # convert nodes to sf
  nodes_sf <- nodes_poly |> 
    arrange(group, x, y) |> 
    group_by(group) |> 
    summarize(
      flowsum = first(flowsum),
      flow = first(flow),
      label = first(label),
      geometry = sf::st_combine(
        sf::st_as_sf(data.frame(x, y), coords = c("x", "y"))
      ) |> 
        sf::st_convex_hull()
    ) |> 
    sf::st_as_sf(crs = crs)
  
  # combine the two sf objects
  nodes_flows <- bind_rows(
    flows_sf |> mutate(type = "flow"),
    nodes_sf |> mutate(type = "node")
  )

  return(nodes_flows)
}

library(yaml)
library(dplyr)
# library(network)
# library(ggnet)
# library(RColorBrewer)
# library(igraph)
library(visNetwork)

concordance_data = read.csv(snakemake@input[["concordance_data"]])
group_annotations = read.csv(snakemake@input[["sample_group_metadata"]], col.names=c("group_id", "pass_fail_type"))
MANIFEST = read.csv(snakemake@input[["annotated_manifest"]])

## Draw Network Plot
#for (group_id in group_annotations$group_id){ 
for (group_id in unique(MANIFEST$group_id)) { 
    print(paste('networkPlot: ', group_id))

    samples_in_group = MANIFEST[MANIFEST$group_id == group_id,] %>% 
      arrange(S_SUBJECTID, SampleTypeCode) %>%
      select(S_SUBJECTID, keep_drop, SampleLabel, SampleTypeCode, SampleAlias, ExpectedGender, ObservedGender, S_SAMPLEID)

    #if (empty(samples_in_group)) {
    #  next
    #}

    starts = c()
    ends = c()
    line_types = c()
    edges <- c()
    c_values <- c()
    for (row1 in 1:nrow(samples_in_group)) {
      id1 = samples_in_group[row1, "SampleAlias"]
      label1 = samples_in_group[row1, "SampleLabel"]
      for (row2 in row1:nrow(samples_in_group)) {
        id2 = samples_in_group[row2, "SampleAlias"]
        label2 = samples_in_group[row2, "SampleLabel"]
        c_obs = concordance_data[concordance_data$ID1==label1 & concordance_data$ID2==label2,] 
        if (nrow(c_obs) == 1) {
            if (c_obs$observed) {
                starts <- append(starts, id1)
                ends <- append(ends, id2)
                line_types <- append(line_types, 'solid')
                c_values <- append(c_values, c_obs$Concord)
            }
        }# else if (nrow(c_obs) > 1) {
        #   c_obs = c_obs[0]
        #   if (c_obs$observed) {
        #       starts <- append(starts, id1)
        #       ends <- append(ends, id2)
        #       line_types <- append(line_types, 'solid')
        #       c_values <- append(c_values, c_obs[0]$Concord)
        #   }
        # }
      }
    }
    print(paste('starts=', starts))
    print(paste('ends=', ends))
    print(paste('c_values=', c_values))
    if (length(starts) > 0) {
      width = 1
      font = "8px courier black"
    } else {
      width = c()
      font = c()
    }

    edges <- data.frame(
      from=starts,
      to=ends,
      title=c_values,
      label=c_values,
      width=width,
      font=font
    )    
    #nodes <- data.frame(
    #  name=samples_in_group$SampleAlias,
    #  carac=samples_in_group$S_SUBJECTID,
    #  frame=samples_in_group$SampleTypeCode
    #)
    subject_ids <- unique(samples_in_group$S_SUBJECTID)
    type_code_shape_map <- data.frame(type=c('d', 'r', 'm'), shape=c("square", "diamond", "triangle"))
    # type_code_icon_map <- data.frame(type=c('d', 'r', 'm'), shape=c("f471", "diamond", "triangle"))
    gender_shape_map <- data.frame(gender=c('M', 'F', 'U', ''), shape=c("square", "dot", "star", "star"))
    palette <- c("#B3E5FC", "#fff176", "#B33F40", "#aed581", "#FFE8E5", "#FFB74D", "#BA68C8", "#A1887F", "#90A4AE", "#E0E0E0", "#FF8A65")

    colors <- data.frame(colors=palette[match(samples_in_group$S_SUBJECTID, subject_ids)])
    print(paste('colors=', colors$colors))

    shadow <- samples_in_group$ObservedGender != samples_in_group$ExpectedGender

    bw <- 2*(samples_in_group$ObservedGender != samples_in_group$ExpectedGender)

    shapes <- gender_shape_map$shape[match(samples_in_group$ObservedGender, gender_shape_map$gender)]

    gender_icon_map <- data.frame(gender=c('M', 'F', 'U'), shape=c("f222", "f221", "f228"))
    gender_icon_map <- data.frame(gender=c('M', 'F', 'U', ''), shape=c("f2a1", "f278", "f100", "f100"))
    icon_codes <- gender_icon_map$shape[match(samples_in_group$ObservedGender, gender_icon_map$gender)]

    nodes <- data.frame(id=samples_in_group$SampleAlias,
                        label=paste0(samples_in_group$SampleAlias, '\n', samples_in_group$S_SAMPLEID, '(', samples_in_group$SampleTypeCode, ')'),
                        group=samples_in_group$S_SUBJECTID,
                        shape=shapes,
                        color=colors$colors,
                        borderWidth=bw,
                        icon.code=icon_codes,
                        icon.color=colors$colors,
                        icon.face="Ionicons",
                        size=17,
                        level=match(samples_in_group$S_SUBJECTID, subject_ids),
                        title=paste0('<b><u>', samples_in_group$S_SUBJECTID,'</u></b><br />',
                                      '<b>Gender=</b>', samples_in_group$ObservedGender, '<br />', 
                                      '<b>SampleType=</b>', samples_in_group$SampleTypeCode, '<br />',
                                      '<b>obs=</b>', samples_in_group$ObservedGender, '<br />',
                                      '<b>exp=</b>', samples_in_group$ExpectedGender, '<br />'
                                    )
                        )
    print(length(subject_ids))
    # png(paste('multiomics/images/group_', group_id, '_network.png'))
    graph = visNetwork(nodes, edges, height="1000px") %>% visLegend() %>% addIonicons #%>% addFontAwesome()
    visSave(plot, paste('multiomics/images/group_', group_id, '_network.html'))
    #dev.off()

}

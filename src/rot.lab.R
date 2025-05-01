rot.lab = function(size, angle){
  if(missing(size)) size = 10
  if(missing(angle)) angle = 40
  theme(axis.text.x = element_text(angle = angle, hjust =1, size = size))
}

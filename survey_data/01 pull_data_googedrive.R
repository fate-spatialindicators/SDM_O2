library(googledrive)

# this should fire up browser window
drive_auth()

target = drive_ls(as_id("https://drive.google.com/drive/u/1/folders/1lD5lrQN01yrmVe_PZa3zIYa8Hm9vQFeX"))

for(i in 1:nrow(target)) {
  # pull in each file as dribble
  drive_download(file=target[i,], 
    path=paste0("survey_data/",target$name[i]), 
    overwrite = TRUE)
}
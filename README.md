# The EDN Suite <img src="https://github.com/angel-bee2018/2021_extended_delta_notation/blob/master/hex_sticker.png" align="right" height="200" />

![code_version](https://img.shields.io/badge/version-pre--publication-ffb3b3) ![known_issues](https://img.shields.io/badge/issues-a%20fucktonne-critical) 

A collection of 4 tools that enable researchers to generate and interpret the Extended Delta notation for alternative splicing. Currently in pre-publication testing.

EDN Automator: Auomatically generates publication-ready Extended Delta Notation to describe human RNA splice events.
EDN Workshop: An extremely powerful tool to visualise the location and consequences of alternative splicing. Plots distances from user exons/junctions to the nearest reference exon and feature vertices. This also includes protein domains from BioMart, Interpro and PTM sites from dbPTM. Users can also upload custom feature tracks.

Available in native R Shiny app (any operating system) or standalone Windows application.

## Standalone app

### Windows

![app_operating_system](https://img.shields.io/badge/standalone-windows_10-a7e1fa)

Link to download the Windows app:

#### Instructions for running Windows standalone app

1. Download and extract the folder called `extended_delta_notation_generator` into the location of your choice. Make sure you have enough storage space as the application is quite large (~6 GB). This is necessary as the program must include all the R libraries and the edited Ensembl GTF in order for R Shiny to work.

2. Open the folder and open `electron-quick-start.exe`. A blank window should appear. Please wait 1 minute as the entire Ensembl GTF needs to be loaded. 

3. Close the blank window and open `electron-quick-start.exe` again. The app should be working now!*

*If anything goes wrong e.g. cannot add extra boxes for exons, then close the window and re-open. This program is not perfect, sorry.

#### Instructions if the program doesn't work and just shows a blank screen

Install electron and set-up the app (this only needs to be done once if running for the first-time):

2. Download npm which should come bundled with node.js here:
    https://nodejs.org/en/download/

3. You need to now install Electron into npm in order for the app to work. To do this, open a Command Prompt in Windows by holding the Windows button and pressing R on the keyboard (Win + R). A box will appear. Type `cmd` into the box. Click "OK". A black box should now appear.

4. Go back and open the folder you just extracted in step 1. Copy the address in the address bar. This is the directory name. In the Command Prompt, type `cd` followed by a space and the directory name, so that it looks like `cd <DIRECTORY>`. Sometimes this can be done quickly by copying the address using (Ctrl + C) and pasting into the Command Prompt by right-clicking the black space in the Command Prompt. This may require a few clicks to work. If it doesn't work, you'll have to manually type the adrress. Sorry.

5. Press enter. The line should now show the address you just entered. If it doesn't work, then type in the letter of the drive that it's located in. So for example, if it's in the C drive, then type `c:` and press enter, before running step 4 again.

6. Into the command prompt, type `npm install electron` and press enter.

7. If it says vulnerabilities were found, follow the instructions to fix these using `npm audit fix`.

8. Now point the terminal to the app folder by typing `cd resources`, enter, followed by `cd app`, enter.

9. Type in `npm start` and press enter. The command prompt will start to print messages beginning with "stderr:", and a blank window should appear. This means that R Shiny is loading.  Please wait a while for it to load the Ensembl annotation file. It should not take more than 10 minutes. When it's finished loading, it should say:

    `stderr:
    Listening on http://127.0.0.1:9191`

10. You can now close the blank window and leave Command Prompt.

11. Return to the open folder and open `electron-quick-start.exe`. That's it!

NOTE: Parallel processing of input GTF files in the Windows standalone app does not work. For parallel processing, please use the native shiny app in R studio instead.

NOTE2: The Delta symbol is not available in the Windows standalone app. The program will therefore output `(Delta)` where `Δ` should be.

Created using R Shiny/Electron.

<img src="https://github.com/rstudio/hex-stickers/blob/master/SVG/shiny.svg" width="100" height="100" /> <img src="https://upload.wikimedia.org/wikipedia/commons/9/91/Electron_Software_Framework_Logo.svg" width="100" height="100" /> 


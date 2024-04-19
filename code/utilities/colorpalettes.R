### general colors used throughout
AMYWT <- "#b2311cff"
AMYKO <- "#e67967ff"

HIPWT <- "#275c62"
HIPKO <- "#47abb5"

P6HIPCTX_XY_WT <- "#4f3f82ff"
P6HIPCTX_XY_KO <- "#806eb8ff"

P6LIVER_WT <- "#a05513ff"
P6LIVER_KO <- "#e18b40ff"

#ESCs
ESC_XY_WT <- "#b68200ff" #(182,130,0)
ESC_XY_KO <- "#f2b929ff" #(242,185,41)

#48 hr epilcs
EpiLC_XY_WT <- "#63344c" #(99,51,76)
EpiLC_XY_KO <- "#b26c93" #(178,108,147)
EpiLC_XX_WT <- "#b84c00ff"
EpiLC_XX_HET <- "#f3873bff"
EpiLC_XX_KO <- "#ffb581ff"

#96 hr epilcs
EpiLC96_XY_WT <- "#134b6e" #(19,75,110)
EpiLC96_XY_KO <- "#3494cb" #(52,147,203)

germcolor <- "#5fa35f"

####### specific palettes #######


#tissue specific genes palettes
tissue_palette <-  c("Adrenal Gland" = "#f5372bff", "Brain" = "#f75c52ff","Forestomach" = "#fa8660ff","Heart" = "#ffa959ff","Kidney" = "#ffd128ff","Liver" = "#edba00ff","Large Intenstine" = "#b7d700ff","Lung" = "#04d700ff","Muscle" = "#39b600ff","Ovary" = "#00c4a5ff","Small Intestine" = "#00cefdff","Spleen" = "#00b9e3ff", "Stomach" = "#5091ffff","Testis" = "#2a7affff","Thymus" = "#8a68ffff","Uterus" = "#db72fbff","Vesicular Gland" = "#ff56b7ff")
tissue_palette2 <-  c("Adrenal Gland" = "#f5372bff", "Brain" = "#f75c52ff","Forestomach" = "#fa8660ff","Heart" = "#ffa959ff","Kidney" = "#ffd128ff","Liver" = "#edba00ff","Large Intenstine" = "#b7d700ff","Lung" = "#04d700ff","Muscle" = "#39b600ff","Ovary" = "#00c4a5ff","Small Intestine" = "#00cefdff","Spleen" = "#00b9e3ff", "Stomach" = "#5091ffff","Testis" = "#2a7affff","Thymus" = "#8a68ffff","Uterus" = "#db72fbff","Vesicular Gland" = "#ff56b7ff", "Non-significant" = "gray38", "Non-specific" = "darkgray")

#all XX and XY EpiLC genos
EpiLCpalette <- c("XYWT" = EpiLC_XY_WT, "XY5cKO" = EpiLC_XY_KO, "XXWT" = EpiLC_XX_WT, "XX5cHET" = EpiLC_XX_HET, "XX5cKO" = EpiLC_XX_KO, "Both" = "deepskyblue2")
#male eplic palette
EpiLC_XY_palette <- c("WT" = EpiLC_XY_WT, "5cKO" = EpiLC_XY_KO, "5CKO" = EpiLC_XY_KO)


#EpiLC with RA
allRA <- c("5CKO0" = ESC_XY_KO, "WT0" = ESC_XY_WT, "5CKO48" = EpiLC_XY_KO, "WT48" = EpiLC_XY_WT, "5CKO96" = EpiLC96_XY_KO, "WT96" = EpiLC96_XY_WT)
RAtimes <- c("0" = ESC_XY_KO, "48" = EpiLC_XY_KO, "96" = EpiLC96_XY_KO)
EpiLCRA_5cKO <- c(EpiLC_XY_KO, EpiLC96_XY_KO)
EpiLCRA_WT <- c("48" = EpiLC_XY_WT, "96" = EpiLC96_XY_WT)
genoRAcolors <- c("WT RA +" = "#d241b2" , "WT RA -" = "#404040", "5CKO RA +" = "#ff85e3" , "5CKO RA -" = "#8f8f8f" )

#palette for EpiLCs and P6 brain
EpiLC_P6Brainpalette <- c("XYWT" = "#283454ff", "XY5cKO" = "#778cbfff", "XXWT" = "#732113ff", "XX5cHET" = "#e36752ff", "XX5cKO" = "#f0ada2ff")

#Sequenced tissues colors
sequenced_tissues_palette <- c("EpiLC" = EpiLC_XY_KO , "P6 Ctx and Hip" = P6HIPCTX_XY_KO, "P6 Liver" = P6LIVER_KO, "Adult Amygdala" = AMYKO, "Adult Hippocampus" = HIPKO)

liverbrain_palette <- c("Both" = germcolor, "Liver" = P6LIVER_KO, "Brain" = P6HIPCTX_XY_KO)
#imprinted genes 
imp_palette <- c("Maternal" = "#fc7cd3" , "Paternal" = "#4aaeff")

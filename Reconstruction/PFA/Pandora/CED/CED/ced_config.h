/********************************************************** 
* ced_config.h, CED config file                           *    
* Hauke Hoelbe, DESY, 2011                                *
* Headerfile to adapt CED before the build.               *
**********************************************************/

#ifndef __CED_CONFIG
#define __CED_CONFIG

#include "ced_cli.h"


/********************************************************** 
* Handling                                                *
**********************************************************/
//enable zoom function by right click and pull
#define ZOOM_RIGHT_CLICK                   0

//time when 2 clicks should be a double click, in 1/1000000 secounds
#define DOUBLE_CLICK_TIME                  300000 

//data layer keys 
#define DATALAYER_SHORTKEY_00       '0'
#define DATALAYER_SHORTKEY_01       '1'
#define DATALAYER_SHORTKEY_02       '2'
#define DATALAYER_SHORTKEY_03       '3'
#define DATALAYER_SHORTKEY_04       '4'
#define DATALAYER_SHORTKEY_05       '5'
#define DATALAYER_SHORTKEY_06       '6'
#define DATALAYER_SHORTKEY_07       '7'
#define DATALAYER_SHORTKEY_08       '8'
#define DATALAYER_SHORTKEY_09       '9'
#define DATALAYER_SHORTKEY_10       ')'
#define DATALAYER_SHORTKEY_11       '!'
#define DATALAYER_SHORTKEY_12       '@'
#define DATALAYER_SHORTKEY_13       '#'
#define DATALAYER_SHORTKEY_14       '$'
#define DATALAYER_SHORTKEY_15       '%'
#define DATALAYER_SHORTKEY_16       '^'
#define DATALAYER_SHORTKEY_17       '&'
#define DATALAYER_SHORTKEY_18       '*'
#define DATALAYER_SHORTKEY_19       '('
#define DATALAYER_SHORTKEY_20       't'
#define DATALAYER_SHORTKEY_21       'y'
#define DATALAYER_SHORTKEY_22       'u'
#define DATALAYER_SHORTKEY_23       'i'
#define DATALAYER_SHORTKEY_24       'o'

//detector layer keys
#define DETECTORLAYER_SHORTKEY_00   'j'
#define DETECTORLAYER_SHORTKEY_01   'k'
#define DETECTORLAYER_SHORTKEY_02   'l'
#define DETECTORLAYER_SHORTKEY_03   ';'
#define DETECTORLAYER_SHORTKEY_04   '\''
#define DETECTORLAYER_SHORTKEY_05   'p'
#define DETECTORLAYER_SHORTKEY_06   '['
#define DETECTORLAYER_SHORTKEY_07   ']'
#define DETECTORLAYER_SHORTKEY_08   '\\'
#define DETECTORLAYER_SHORTKEY_09   'T'
#define DETECTORLAYER_SHORTKEY_10   'Y'
#define DETECTORLAYER_SHORTKEY_11   'U'
#define DETECTORLAYER_SHORTKEY_12   'I'
#define DETECTORLAYER_SHORTKEY_13   'O'
#define DETECTORLAYER_SHORTKEY_14   'P'
#define DETECTORLAYER_SHORTKEY_15   '{'
#define DETECTORLAYER_SHORTKEY_16   '}'
#define DETECTORLAYER_SHORTKEY_17   '|'
#define DETECTORLAYER_SHORTKEY_18   'a'
#define DETECTORLAYER_SHORTKEY_19   'e'


/********************************************************** 
* Colors and appearance                                   *
**********************************************************/
//size of the boarder line in filled (new view) of detector components.
#define CED_GEOTUBE_LINE_WIDTH              0.3

//maximal transparency of boarder lines 
//#define CED_GEOTUBE_LINE_MAX_TRANS          0.2  
#define CED_GEOTUBE_LINE_MAX_TRANS          1.0  


//names and values of color apairs in popup menu
#define CED_BGCOLOR_OPTION1_NAME            "Gainsboro" 
#define CED_BGCOLOR_OPTION1_COLORCODE       0.862745,0.862745,0.862745,0 

#define CED_BGCOLOR_OPTION2_NAME            "Lightgrey" 
#define CED_BGCOLOR_OPTION2_COLORCODE       0.827451,0.827451,0.827451,0

#define CED_BGCOLOR_OPTION3_NAME            "Darkgray" 
#define CED_BGCOLOR_OPTION3_COLORCODE       0.662745,0.662745,0.662745,0

#define CED_BGCOLOR_OPTION4_NAME            "Gray" 
#define CED_BGCOLOR_OPTION4_COLORCODE       0.501961,0.501961,0.501961,0

#define CED_BGCOLOR_OPTION5_NAME            "Silver" 
#define CED_BGCOLOR_OPTION5_COLORCODE       0.7529,0.7529,0.7529,0

#define CED_BGCOLOR_OPTION6_NAME            "Dimgray" 
#define CED_BGCOLOR_OPTION6_COLORCODE       0.4118,0.4118,0.4118,0

#define CED_BGCOLOR_OPTION7_NAME            "Lightsteelblue" 
#define CED_BGCOLOR_OPTION7_COLORCODE       0.6902,0.7686 ,0.8706,0

#define CED_BGCOLOR_OPTION8_NAME            "Steelblue" 
#define CED_BGCOLOR_OPTION8_COLORCODE       0.2745,0.5098,0.70588,0

#define CED_BGCOLOR_OPTION9_NAME            "Seagreen" 
#define CED_BGCOLOR_OPTION9_COLORCODE       0.18039,0.54509,0.34117,0

#define CED_BGCOLOR_OPTION10_NAME           "Orange" 
#define CED_BGCOLOR_OPTION10_COLORCODE      1,0.647,0,0

#define CED_BGCOLOR_OPTION11_NAME           "Yellow" 
#define CED_BGCOLOR_OPTION11_COLORCODE      1,1,0,0

#define CED_BGCOLOR_OPTION12_NAME           "Violett" 
#define CED_BGCOLOR_OPTION12_COLORCODE      0.9333,0.5098,0.9333,0

#define CED_BGCOLOR_OPTION13_NAME           "Black" 
#define CED_BGCOLOR_OPTION13_COLORCODE      0,0,0,0

#define CED_BGCOLOR_OPTION14_NAME           "Blue" 
#define CED_BGCOLOR_OPTION14_COLORCODE      0,0.2,0.4,0

#define CED_BGCOLOR_OPTION15_NAME           "White" 
#define CED_BGCOLOR_OPTION15_COLORCODE      1,1,1,0




//Color of xyz axes 
#define AXES_COLOR                          0.2,0.2,0.8

//Width of xyz axes 
#define AXES_LINE_SIZE                      0.5


//Help frame: Frame fill color, and transp
#define HELP_FRAME_FILL_COLOR               0.5,1,1,0.8

//Help frame: Frame boarder color and transp
#define HELP_FRAME_BOARDER_COLOR            0.1,0.8,1.0,0.8 

//Help frame: Frame boarder line width
#define HELP_FRAME_BOARDER_LINE_SIZE        3. 

//Help frame: Text color, and transp
#define HELP_FRAME_TEXT_COLOR               0.0,0.0,0.0 




/********************************************************** 
* Layers                                                  *
**********************************************************/
//number of total number of layers 
#define CED_MAX_LAYER                       100

//number of layers shown in popup menu
#define NUMBER_POPUP_LAYER                  20

//number of layers reserved for data
#define NUMBER_DATA_LAYER                   25

//number of layers reserved for detector components 
#define NUMBER_DETECTOR_LAYER               20

//layer description text: maximal number of chars for one entry
#define CED_MAX_LAYER_CHAR                  400

/********************************************************** 
* Graphics                                                *
**********************************************************/

//Camera field of view, in degree 
#define CAMERA_FIELD_OF_VIEW                45

//Camera min distance (hint: min and max should be close together)
#define CAMERA_MIN_DISTANCE                 100

//Camera max distance (hint: min and max should be close together)
#define CAMERA_MAX_DISTANCE                 50000.0*mm.sf+50000/mm.sf 

//Where the camera stands
#define CAMERA_POSITION                     0,0,2000



//Fisheye alpha factor
#define FISHEYE_ALPHA                       1e-3
#define FISHEYE_ZOOM                        8.

/********************************************************** 
* Debug                                                   *
**********************************************************/

//show pickable points:
//#define DEBUG_PICKING 1




#endif 

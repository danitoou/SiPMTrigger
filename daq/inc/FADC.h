#ifndef _FADC_H_
#define _FADC_H_

int FADCTriggerConfig();

int FADCSetChannelOffsets();


int FADCBoardInit();

int FADCIOSetup();

int FADCDigitizationSetup();

int FADCBoardExit();

int FADCExit();

int FADCStartAcquisition();
int FADCStopAcquisition();

int FADCSetChannelTriggers();
int FADCCleanErrors();
int FADCCleanBuffers();
int FADCChannelsCallibrate();

int FADCPostConfig();


#endif

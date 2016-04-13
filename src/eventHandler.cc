#include "eventHandler.h"






eventHandler::~eventHandler() {


}

void eventHandler::print() {
  std::cout << "Print me!" << std::endl;

}


void eventHandler::eventSlot(Int_t event, Int_t x, Int_t y, TObject *selected) {

  switch (event) {
  case 1: { 
    Double_t px;
    Double_t py;
    tc.AbsPixeltoXY(x, y, px, py);
    std::cout << "Setting Threshold to: " << py << std::endl; 
    _threshold = py; 

    break; 
  }

  default:
    break;
    //    std::cout << "Event ID:" << event << " X:" << x << " Y:" << y << std::endl;
  }

  
}



  

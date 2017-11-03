//
// Created by brad on 3/9/16.
//

#include "ZREval.h"

int main(int argc, char* argv[]) {

  int lang=0; //0 is english 1 is tsonga
  std::string path = "/home/brad/data/english_wavs";
  //std::string path = "/home/brad/data/tsonga";

  if (argc > 1) {
    path = std::string(argv[1]);
  }

  if (argc > 2) {
    lang = atoi(argv[2]);
  }

  ZREval evaluation(lang);

  evaluation.runEvalFromCfg(path);

  return 0;
}

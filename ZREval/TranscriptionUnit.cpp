//
// Created by brad on 3/19/16.
//
#include "DataTypes.h"

int minTUSize = 480; //30 ms

void TranscriptionUnit::insertClass(Class newMember) {
  //search through subsequences for case where ss.start >= new Member.start && ss.end <= new member.end
  //might want to check abour the effect averaging has on this

  //for(segment s: subSequences){
  int indexToBeRplaced = 0;
  for (int i = 0; i < subSequences.size(); i++) {
    //check if new member ris within bounds
    if (subSequences[i].start >= newMember.start && subSequences[i].end <= newMember.end) {

      int origStart = subSequences[i].start;
      int origEnd = subSequences[i].end;

      bool completelyCovered = true;

      //check if there is enough room before inserted segment
      if (newMember.start - origStart >= minTUSize) {
        subSequences[i].end = newMember.start - 1;
        completelyCovered = false;
      }

      //if there is no remainder before insert newMember in place of current segment
      if (completelyCovered) {
        subSequences.insert(subSequences.begin() + i, newMember);
      }
      else { //if there is a remainder before, insert it after
        subSequences.insert(subSequences.begin() + i + 1, newMember);
      }

      //check if there is space at the end
      if (origEnd - newMember.end >= minTUSize) {
        //if we already generated a remainder for the first part we need to make a new class
        if (!completelyCovered) {
          Class endChunk(newMember.end + 1, origEnd, false);
          subSequences.insert(subSequences.begin()+i+2, endChunk);
        }
        else { //if we didn't find one before it we don't need to make a new one
          subSequences[i].start = newMember.end + 1;
          completelyCovered = false;
        }
      }

      //if the previously existing segment is completely covered by the new member, delete it
      if (completelyCovered) {
        //it's i+1 because we already inserted newMember
        subSequences.erase(subSequences.begin() + i + 1);
      }
      break; //we dont care about there rest of the subsequences, because they won't overlap. hopefully
    }
  }
}
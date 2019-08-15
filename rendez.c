#include <stdio.h>
#include <stdlib.h>

/*
** CONSTS **
*/
#define ZMAX 86

const char EL[ZMAX][4] = {
  "H" ,                                                                                                                                                                                     "He",
  "Li", "Be", "B" ,                                                                                                                                                 "C" , "N" , "O" , "F" , "Ne",
  "Na", "Mg", "Al",                                                                                                                                                 "Si", "P" , "S" , "Cl", "Ar",
  "K" , "Ca", "Sc",                                                                                     "Ti", "V" , "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
  "Rb", "Sr", "Y" ,                                                                                     "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I" , "Xe",
  "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W" , "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn"
};

int getElId(char *str, int *i)
{
  int id = 0;
  int start = *i;
  while(EL[id][0] != str[*i])
    id++;
  (*i)++;
  while(str[*i] >= 'a' && str[*i] <= 'z')
  {
    while(EL[id][(*i) - start] != str[*i])
      id++;
    (*i)++;
  }
  return id;
}

/*
** STRUCTS AND PARSING **
*/
struct Field
{
  int element;
  int count;
};

void readField(char *str, int *i, struct Field *field)
{
  field->count = 0;
  field->element = getElId(str, i);
  while(str[*i] >= '0' && str[*i] <= '9')
  {
    field->count *= 10;
    field->count += str[*i] - '0';
    (*i)++;
  }
  if(field->count == 0)
    field->count = 1;
}

struct Molecule
{
  int fieldCount;
  struct Field fields[16];
  int charge;
  int count;
};

char readMolecule(char *str, int *i, struct Molecule *molecule)
{
  molecule->fieldCount = 0;
  molecule->charge = 0;
  molecule->count = 1;
  while(str[*i] > 'A' && str[*i] < 'Z')
  {
    readField(str, i, molecule->fields + molecule->fieldCount);
    molecule->fieldCount++;
  }
  if(str[*i] == '^')
  {
    molecule->charge = str[(*i) + 3] - '0';
    if(str[(*i) + 2] == '-')
      molecule->charge *= -1;
    *i += 5;
  }
  (*i)++;
  return str[(*i) - 1];
}

/*
** DATA **
*/
int leftCount = 0;
struct Molecule left[16];
int rightCount = 0;
struct Molecule right[16];

/*
** SOLVING **
*/
int fieldCount = 0;
struct Field fields[16];
int rightCounts[16];

void printResult()
{
    printf("Success, the coefficients are:\n");
    for(int j2 = 0; j2 < leftCount; j2++)
    {
      printf("%d ", left[j2].count);
    }
    printf("\n");
    for(int j2 = 0; j2 < rightCount; j2++)
    {
      printf("%d ", right[j2].count);
    }
    printf("\n");
}

int genRight(int j)
{
  int bonus = 0;
  int over = 0;
  while(1)
  {
    int success = 1;
    for(int l = 0; l < fieldCount; l++)
    {
      if(rightCounts[l] != fields[l].count)
      {
        success = 0;
        break;
      }
    }
    if(success)
    {
      printResult();
      return 1;
    }
    if(j < rightCount - 1)
      if(genRight(j + 1))
        return 1;
    for(int k = 0; k < right[j].fieldCount; k++)
      for(int l = 0; l < fieldCount; l++)
        if(fields[l].element == right[j].fields[k].element)
        {
          rightCounts[l] += right[j].fields[k].count;
          if(rightCounts[l] > fields[l].count)
            over = 1;
          break;
        }
    right[j].count++;
    bonus++;
    if(over)
    {
      for(int k = 0; k < right[j].fieldCount; k++)
        for(int l = 0; l < fieldCount; l++)
          if(fields[l].element == right[j].fields[k].element)
          {
            rightCounts[l] -= right[j].fields[k].count * bonus;
            break;
          }
      right[j].count -= bonus;
      return 0;
    }
  }
}

int match()
{
  fieldCount = 0;
  for(int j = 0; j < leftCount; j++)
  {
    for(int k = 0; k < left[j].fieldCount; k++)
    {
      int newElement = 1;
      for(int l = 0; l < fieldCount; l++)
      {
        if(fields[l].element == left[j].fields[k].element)
        {
          fields[l].count += left[j].fields[k].count * left[j].count;
          newElement = 0;
          break;
        }
      }
      if(newElement)
      {
        fields[fieldCount].element = left[j].fields[k].element;
        fields[fieldCount].count = left[j].fields[k].count * left[j].count;
        fieldCount++;
      }
    }
  }
  for(int l = 0; l < fieldCount; l++)
    rightCounts[l] = 0;
  for(int j = 0; j < rightCount; j++)
    for(int k = 0; k < right[j].fieldCount; k++)
      for(int l = 0; l < fieldCount; l++)
        if(fields[l].element == right[j].fields[k].element)
        {
          rightCounts[l] += right[j].fields[k].count;
          if(rightCounts[l] > fields[l].count)
            return 0;
          break;
        }
  return genRight(0);
}

int confLeft(int bonus, int step, int start)
{
  if(step >= bonus)
    return match();
  for(int j = start; j < leftCount; j++)
  {
    left[j].count++;
    if(confLeft(bonus, step + 1, j))
      return 1;
    left[j].count--;
  }
  return 0;
}

/*
** MAIN **
*/
int main(int argc, char **argv)
{
  if(argc != 2)
  {
    fprintf(stderr, "Invalid arguments.\n");
    return 1;
  }
  int i = 0;
  while(1)
  {
    char n = readMolecule(argv[1], &i, left + leftCount);
    leftCount++;
    if(n == '=')
      break;
    if(n != '+')
    {
      fprintf(stderr, "Invalid equation.\n");
      return 1;
    }
  }
  while(1)
  {
    char n = readMolecule(argv[1], &i, right + rightCount);
    rightCount++;
    if(n == 0)
      break;
    if(n != '+')
    {
      fprintf(stderr, "Invalid equation.\n");
      return 1;
    }
  }
  for(int bonus = 0; 1; bonus++)
  {
    printf("Trying bonus %d\n", bonus);
    if(confLeft(bonus, 0, 0))
      return 0;
  }
}


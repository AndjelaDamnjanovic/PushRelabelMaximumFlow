#include "Headers/graphwindow.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    GraphWindow w;
    w.show();
    return a.exec();
}

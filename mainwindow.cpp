#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "mainwidget.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    QWidget *mainwidget = new MainWidget(this);
    this->setCentralWidget(mainwidget);
    this->setMinimumSize(400, 400);
}

MainWindow::~MainWindow()
{
    delete ui;
}

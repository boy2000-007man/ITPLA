#ifndef MAINWIDGET_H
#define MAINWIDGET_H

#include <QWidget>

class MainWidget : public QWidget
{
    Q_OBJECT
public:
    explicit MainWidget(QWidget *parent = 0);
private:
    void paintEvent(QPaintEvent *);
signals:

public slots:

};

#endif // MAINWIDGET_H

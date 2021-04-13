# goodbye server

observeEvent(input$exit, {
  js$closewindow();
  stopApp(7)
})
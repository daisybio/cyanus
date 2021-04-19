# goodbye server

observeEvent(input$exit, {
  js$closewindow();
  stopApp("Quit App on Goodbye...")
})
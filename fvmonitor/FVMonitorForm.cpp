#include "FVMonitorForm.h"

using namespace System;
using namespace System::Windows::Forms;


[STAThread]
void Main(array<String^>^ args)
{
	Application::EnableVisualStyles();
	Application::SetCompatibleTextRenderingDefault(false);

	fvmonitor::FVMonitorForm form;
	Application::Run(%form);
}
